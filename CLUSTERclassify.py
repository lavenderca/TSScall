#!/usr/bin/env python

# CREATED BY CHRISTOPHER LAVENDER
# BASED ON WORK BY ADAM BURKHOLDER
# INTEGRATIVE BIOINFORMATICS, NIEHS
# WORKING OBJECT ORIENTED VERSION

import os
import argparse
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
            print(line.strip().split('\t')[0:16])
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
                    if annotation[key][0] is not None:
                        annotation_value = annotation[key][0]
                    else:
                        annotation_value = annotation['gene_id']
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


class ClusterClassify(object):

    def __init__(self, **kwargs):

        self.input_detail_file = kwargs['input_detail_file']
        self.annotation_file = kwargs['annotation_file']
        self.output_file = kwargs['output_file']
        self.id_field = kwargs['identification_field']
        self.proximity_threshold = kwargs['proximity_threshold']

        assert os.path.exists(self.input_detail_file)
        assert os.path.exists(self.annotation_file)

        self.execute()

    def printOutput(self, clusters, output_file):
        cluster_lines = dict()
        with open(output_file, 'w') as OUTPUT:
            OUTPUT.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
                         .format(
                            'Cluster ID',
                            'TSSs in cluster',
                            'Chromosome',
                            'Cluster start',
                            'Cluster end',
                            'Overlapping',
                            'Proximal (within ' +
                            str(self.proximity_threshold) + ')',
                            'Representative TSS',
                            'Representative TSS position',
                            'Representative TSS strand',
                            'Closest',
                            'Distance to closest',
                         ))
            for chrom in clusters:
                for cluster in sorted(clusters[chrom],
                                      key=itemgetter('cluster_id')):
                    overlapping = ';'.join(cluster['overlapping'])
                    proximal = ';'.join(cluster['proximal'])
                    line = (
                        '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
                        .format(
                                cluster['cluster_id'],
                                cluster['cluster_members'],
                                cluster['chromosome'],
                                str(cluster['start']),
                                str(cluster['end']),
                                overlapping,
                                proximal,
                                cluster['representative_tss_id'],
                                cluster['representative_tss_position'],
                                cluster['representative_tss_strand'],
                                cluster['closest_id'],
                                cluster['closest_dist'],
                        ))
                    cluster_lines[cluster['cluster_id']] = line
            for cluster_id in sorted(cluster_lines):
                OUTPUT.write(cluster_lines[cluster_id])

    def execute(self):
        annotation = readInReferenceAnnotation(self.annotation_file)[0]
        annotation_ranges = makeRangesFromAnnotation(annotation)
        cluster_ranges = readInClusters(self.input_detail_file)

        findOverlaps(cluster_ranges, annotation_ranges, self.id_field,
                     self.proximity_threshold)

        self.printOutput(cluster_ranges, self.output_file)

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
