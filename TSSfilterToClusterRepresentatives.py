#!/usr/bin/env python

# CREATED BY CHRISTOPHER LAVENDER
# INTEGRATIVE BIOINFORMATICS, NIEHS

import argparse


def printBedEntry(entry, OUTPUT):
    OUTPUT.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        entry['chromosome'],
        str(int(entry['position'])-1),
        entry['position'],
        entry['tss_id'],
        '0',
        entry['strand'],
        ))


def selectClusterRepresentatives(input_detail, output_bed, output_detail):
    clusters = dict()
    detail_entries = dict()

    if output_detail:
        DETAIL_OUTPUT = open(output_detail, 'w')

    with open(input_detail) as f, open(output_bed, 'w') as BED_OUTPUT:
        header = next(f)
        if output_detail:
            DETAIL_OUTPUT.write(header)
        for line in f:
            [
                tss_id,
                tss_type,
                trs,
                genes,
                strand,
                chromosome,
                position,
                reads,
                divergent_flag,
                divergent_partner,
                divergent_distance,
                convergent_flag,
                convergent_partner,
                convergent_distance,
                tss_cluster,
                cluster_members
            ] = line.strip().split('\t')[0:16]

            detail_entries[tss_id] = line

            if tss_cluster != "NA":
                cluster_entry = {
                        'tss_id': tss_id,
                        'tss_type': tss_type,
                        'strand': strand,
                        'chromosome': chromosome,
                        'position': position,
                        'reads': reads,
                }

                if tss_cluster in clusters:
                    clusters[tss_cluster].append(cluster_entry)
                else:
                    clusters[tss_cluster] = [cluster_entry]

        for cluster in sorted(clusters):
            max_reads = 0
            max_read_entry = None
            for entry in clusters[cluster]:
                if entry['tss_type'] == 'called from reference window':
                    max_reads = float('Inf')
                    max_read_entry = entry
                elif int(entry['reads']) > max_reads:
                    max_reads = int(entry['reads'])
                    max_read_entry = entry
                elif int(entry['reads']) == max_reads:
                    if entry['tss_id'] < max_read_entry['tss_id']:
                        max_read_entry = entry
            printBedEntry(max_read_entry, BED_OUTPUT)
            if output_detail:
                DETAIL_OUTPUT.write(detail_entries[max_read_entry['tss_id']])

    if output_detail:
        DETAIL_OUTPUT.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_detail', type=str, default=None,
                        help='output detail file with cluster representatives')
    parser.add_argument('input_detail_file', type=str,
                        help='detail file from TSScall.py')
    parser.add_argument('output_bed', type=str,
                        help='output BED file with cluster representatives')
    args = parser.parse_args()

    selectClusterRepresentatives(args.input_detail_file, args.output_bed,
                                 args.output_detail)
