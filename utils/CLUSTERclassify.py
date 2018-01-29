#!/usr/bin/env python

# CREATED BY CHRISTOPHER LAVENDER
# BASED ON WORK BY ADAM BURKHOLDER
# INTEGRATIVE BIOINFORMATICS, NIEHS
# WORKING OBJECT ORIENTED VERSION

import os
import argparse
import csv
from operator import itemgetter
from TSScall import readInReferenceAnnotation
from collections import defaultdict


def makeRangesFromAnnotation(annotation):
    ranges = defaultdict(list)
    for entry in annotation:
        ranges[annotation[entry]['chromosome']].append({
            'transcript_id': entry,
            'chromosome': annotation[entry]['chromosome'],
            'start': annotation[entry]['tr_start'],
            'end': annotation[entry]['tr_end'],
            'strand': annotation[entry]['strand'],
            'gene_id': annotation[entry]['gene_id'],
            'tss': annotation[entry]['tss'],
        })
        for field in annotation[entry]['gtf_fields']:
            ranges[annotation[entry]['chromosome']][-1][field] = \
                annotation[entry]['gtf_fields'][field]
    return ranges


def readInClusters(input_detail_file):
    clusters = dict()
    cluster_ranges = defaultdict(list)

    with open(input_detail_file) as f:
        next(f)
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

            if tss_cluster != "NA":
                cluster_entry = {
                        'tss_id': tss_id,
                        'tss_type': tss_type,
                        'strand': strand,
                        'chromosome': chromosome,
                        'position': int(position),
                        'reads': int(reads),
                        'cluster_members': cluster_members,
                }

                if tss_cluster in clusters:
                    clusters[tss_cluster].append(cluster_entry)
                else:
                    clusters[tss_cluster] = [cluster_entry]

    for cluster in clusters:
        max_reads = 0
        max_read_entry = None

        cluster_start = float('Inf')
        cluster_end = float('-Inf')

        for entry in clusters[cluster]:
            if entry['tss_type'] == 'called from reference window':
                max_reads = float('Inf')
                max_read_entry = entry
            elif entry['reads'] > max_reads:
                max_reads = entry['reads']
                max_read_entry = entry
            elif entry['reads'] == max_reads:
                if entry['tss_id'] < max_read_entry['tss_id']:
                    max_read_entry = entry
            if cluster_start > entry['position']:
                cluster_start = entry['position']
            if cluster_end < entry['position']:
                cluster_end = entry['position']

        cluster_ranges[max_read_entry['chromosome']].append({
            'cluster_id': cluster,
            'representative_tss_id': max_read_entry['tss_id'],
            'representative_tss_position': max_read_entry['position'],
            'representative_tss_strand': max_read_entry['strand'],
            'cluster_members': max_read_entry['cluster_members'],
            'start': cluster_start,
            'end': cluster_end,
            'chromosome': max_read_entry['chromosome'],
        })

    return cluster_ranges


def findOverlaps(clusters, annotations, key, window):
    def checkOverlap(entry_1, entry_2, threshold):
        if entry_1['chromosome'] == entry_2['chromosome']:
            if entry_1['start'] - threshold <= entry_2['start'] and\
                    entry_1['end'] + threshold >= entry_2['start']:
                return True
            elif entry_1['start'] - threshold <= entry_2['end'] and\
                    entry_1['end'] + threshold >= entry_2['end']:
                return True
            elif entry_1['start'] - threshold >= entry_2['start'] and\
                    entry_1['end'] + threshold <= entry_2['end']:
                return True
        return False

    for chrom in clusters:
        for cluster in clusters[chrom]:
            overlaps = []
            proximal = []
            closest_dist = float('Inf')
            closest_tss = None
            closest_value = None

            for annotation in annotations[chrom]:
                if cluster['chromosome'] == annotation['chromosome']:
                    if annotation[key][0] is not None and key != 'gene_id' \
                            and key != 'transcript_id':
                        annotation_value = annotation[key][0]
                    else:
                        annotation_value = annotation[key]
                    if checkOverlap(cluster, annotation, 0):
                        if annotation_value not in overlaps:
                            overlaps.append(annotation_value)
                    elif checkOverlap(cluster, annotation, window):
                        if annotation_value not in overlaps and\
                                annotation_value not in proximal:
                            proximal.append(annotation_value)

                    tss_distance = abs(cluster['representative_tss_position'] -
                                       annotation['tss'])
                    if closest_tss is None:
                        closest_tss = annotation
                        closest_value = annotation_value
                        closest_dist = tss_distance
                    elif tss_distance < closest_dist:
                        closest_tss = annotation
                        closest_value = annotation_value
                        closest_dist = tss_distance
                    elif tss_distance == closest_dist:
                        # tie-breakers: (1) upstream?, (2) plus strand?,
                        # (3) name
                        if annotation['tss'] < closest_tss['tss']:
                            closest_tss = annotation
                            closest_value = annotation_value
                        elif annotation['tss'] == closest_tss['tss']:
                            if annotation['strand'] >\
                                    closest_tss['strand']:
                                closest_tss = annotation
                                closest_value = annotation_value
                            elif annotation['strand'] ==\
                                    closest_tss['strand']:
                                if annotation_value < closest_value:
                                    closest_value = annotation_value

            cluster.update({
                'overlapping': overlaps,
                'proximal': proximal,
            })

            if closest_tss:
                cluster.update({
                    'closest_id': closest_value,
                    'closest_dist': closest_dist,
                })
            else:
                cluster.update({
                    'closest_id': 'NA',
                    'closest_dist': 'NA',
                })


def findClosestActiveTSS(input_detail_file):
    closest_tss = dict()

    with open(input_detail_file) as f:
        tss_list = csv.DictReader(f, skipinitialspace=True, delimiter='\t')
        sorted_list = sorted(tss_list, key=lambda k: (
            k['Chromosome'],
            k['Position']
            ))

        for i, tss in enumerate(sorted_list):

            upstream_tss = None
            downstream_tss = None

            index = i
            while index >= 0 and upstream_tss is None:
                if sorted_list[index]['Type'] == \
                        'called from reference window' \
                        and tss['Chromosome'] == \
                        sorted_list[index]['Chromosome']:
                    upstream_tss = sorted_list[index]
                index -= 1

            index = i
            while index < len(sorted_list) and downstream_tss is None:
                if sorted_list[index]['Type'] == \
                        'called from reference window' \
                        and tss['Chromosome'] == \
                        sorted_list[index]['Chromosome']:
                    downstream_tss = sorted_list[index]
                index += 1

            if upstream_tss and downstream_tss:
                upstream_distance = abs(
                    int(tss['Position']) - int(upstream_tss['Position']))
                downstream_distance = abs(
                    int(tss['Position']) - int(downstream_tss['Position']))
                if upstream_distance < downstream_distance:
                    closest = upstream_tss
                elif downstream_distance < upstream_distance:
                    closest = downstream_tss
                elif upstream_tss == downstream_tss:
                    closest = upstream_tss
            elif upstream_tss:
                closest = upstream_tss
            elif downstream_tss:
                closest = downstream_tss
            else:
                closest = None

            closest_tss[tss['TSS ID']] = closest
    return closest_tss


class ClusterClassify(object):

    def __init__(self, **kwargs):

        self.input_detail_file = kwargs['input_detail_file']
        self.annotation_file = kwargs['annotation_file']
        self.output_file = kwargs['output_file']
        self.id_field = kwargs['identification_field']
        self.proximity_threshold = kwargs['proximity_threshold']

        assert os.path.exists(self.input_detail_file)
        assert os.path.exists(self.annotation_file)

        self.non_attribute_fields = [
            'TSS ID',
            'Type',
            'Transcripts',
            'Gene ID',
            'Strand',
            'Chromosome',
            'Position',
            'Reads',
            'Divergent?',
            'Divergent partner',
            'Divergent distance',
            'Convergent?',
            'Convergent partner',
            'Convergent distance',
            'TSS cluster',
            'TSSs in associated cluster',
        ]
        self.execute()

    def printOutput(self, clusters, closest, output_file):

        attribute_fields = []
        for cluster_id in closest:
            if closest[cluster_id] is not None:
                for key in closest[cluster_id]:
                    if key not in self.non_attribute_fields and \
                            key not in attribute_fields:
                        attribute_fields.append(key)

        cluster_lines = dict()
        with open(output_file, 'w') as OUTPUT:
            # OUTPUT.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
            #              .format(
            header_fields = [
                'Cluster ID',
                'TSSs in cluster',
                'Chromosome',
                'Cluster start',
                'Cluster end',
                'Representative TSS',
                'Representative TSS position',
                'Representative TSS strand',
                'Overlapping genes',
                'Proximal genes (within ' +
                str(self.proximity_threshold) + ')',
                'Closest gene',
                'Distance to closest gene',
                'Closest active TSS ID',
                'Distance to closest active TSS',
                'Closest active TSS transcript ID',
                'Closest active TSS gene ID',
                'Closest active TSS chromosome',
                'Closest active TSS position',
                'Closest active TSS strand',
            ]
            for field in attribute_fields:
                header_fields.append('Closest active TSS ' + field.lower())
            for i, field in enumerate(header_fields):
                if i == 0:
                    OUTPUT.write(field)
                else:
                    OUTPUT.write('\t' + field)
            OUTPUT.write('\n')

            for chrom in clusters:
                for cluster in sorted(clusters[chrom],
                                      key=itemgetter('cluster_id')):
                    overlapping = ';'.join(cluster['overlapping'])
                    proximal = ';'.join(cluster['proximal'])
                    line = (
                        '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'
                        .format(
                                cluster['cluster_id'],
                                cluster['cluster_members'],
                                cluster['chromosome'],
                                str(cluster['start']),
                                str(cluster['end']),
                                cluster['representative_tss_id'],
                                cluster['representative_tss_position'],
                                cluster['representative_tss_strand'],
                                overlapping,
                                proximal,
                                cluster['closest_id'],
                                cluster['closest_dist'],
                        ))

                    tss_id = cluster['representative_tss_id']
                    if closest[tss_id] is not None:
                        line += ('\t{}' * 7).format(
                            closest[tss_id]
                            ['TSS ID'],
                            str(abs(int(closest[tss_id]['Position']) -
                                int(cluster['representative_tss_position']))),
                            closest[tss_id]['Transcripts'],
                            closest[tss_id]['Gene ID'],
                            closest[tss_id]['Chromosome'],
                            closest[tss_id]['Position'],
                            closest[tss_id]['Strand'],
                        )
                        for field in attribute_fields:
                            if field in closest[tss_id]:
                                line += '\t' + closest[tss_id][field]
                            else:
                                line += '\tNA'
                    else:
                        for i in range(6 + len(attribute_fields)):
                            line += '\tNA'

                    line += '\n'

                    cluster_lines[cluster['cluster_id']] = line
            for cluster_id in sorted(cluster_lines):
                OUTPUT.write(cluster_lines[cluster_id])

    def execute(self):
        annotation = readInReferenceAnnotation(self.annotation_file)[0]
        annotation_ranges = makeRangesFromAnnotation(annotation)
        cluster_ranges = readInClusters(self.input_detail_file)

        closest_TSSs = findClosestActiveTSS(self.input_detail_file)
        findOverlaps(cluster_ranges, annotation_ranges, self.id_field,
                     self.proximity_threshold)

        self.printOutput(cluster_ranges, closest_TSSs, self.output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--proximity_threshold', '-p', default=10000, type=int,
                        help='distance threshold used to determine proximity\
                         (Default: 10000)')
    parser.add_argument('--identification_field', '-i', default='gene_id',
                        type=str, help='field in GTF used to differentiate and\
                        identify genes (Default: "gene_id")')
    parser.add_argument('input_detail_file', type=str,
                        help='input detail file from TSScall.py')
    parser.add_argument('annotation_file', type=str,
                        help='annotation in GTF format')
    parser.add_argument('output_file', type=str,
                        help='output tab-delimited cluster information file')
    args = parser.parse_args()

ClusterClassify(**vars(args))
