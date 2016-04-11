#!/usr/bin/env python

## CREATED BY CHRISTOPHER LAVENDER
## BASED ON WORK BY ADAM BURKHOLDER
## INTEGRATIVE BIOINFORMATICS, NIEHS
## WORKING OBJECT ORIENTED VERSION

import sys
import os
import math
import argparse
import operator

class TSSCalling(object):

    def __init__(self, **kwargs):

        self.forward_bedgraph = kwargs['forward_bedgraph']
        self.reverse_bedgraph = kwargs['reverse_bedgraph']
        self.chrom_sizes = kwargs['chrom_sizes']
        self.annotation_file = kwargs['annotation_file']
        self.output_bed = kwargs['output_bed']

        assert os.path.exists(self.forward_bedgraph)
        assert os.path.exists(self.reverse_bedgraph)
        assert os.path.exists(self.chrom_sizes)
        if self.annotation_file:
            assert os.path.exists(self.annotation_file)

        self.fdr_threshold = kwargs['fdr']
        self.false_positives = kwargs['false_positives']
        self.nutss_filter_size = kwargs['nutss_filter_size']
        self.nutss_search_window = kwargs['nutss_search_window']
        self.bidirectional_threshold = kwargs['bidirectional_threshold']
        self.cluster_threshold = kwargs['cluster_threshold']
        self.detail_file = kwargs['detail_file']

        self.tss_list = []
        self.reference_annotation = None
        self.annotated_tss_count = 0
        self.unannotated_tss_count = 0
        self.tss_cluster_count = 0

        self.execute()

    ## STRAND_STATUS IS USED TO DETERMINE IF STRAND IS USED IN SORT
    def sortList(self, input_list, strand_status):
        if strand_status == 'sort_by_strand':
            return sorted(input_list, key=lambda k: (
                k['strand'],
                k['chromosome'],
                k['start']
                ))
        elif strand_status == 'ignore_strand':
            return sorted(input_list, key=lambda k: (
                k['chromosome'],
                k['start']
                ))

    ## ENTRY 1 IS LESS THAN ENTRY 2?
    def isLessThan(self, entry_1, entry_2):
        for val in ['strand', 'chromosome', 'start']:
            if entry_1[val] < entry_2[val]:
                return True
            elif entry_1[val] > entry_2[val]:
                return False
        return False
    ## ENTRY 1 IS WITHIN ENTRY 2?
    def isWithin(self, entry_1, entry_2):
        if entry_1['strand'] == entry_2['strand'] and entry_1['chromosome'] == entry_2['chromosome']:
            if entry_1['start'] >= entry_2['start'] and entry_1['end'] <= entry_2['end']:
                return True
        return False

    def getID(self, base_name, count):
        max_entries = 999999
        feature_name = base_name + '_'
        for i in range(len(str(count)), len(str(max_entries))):
            feature_name += '0'
        feature_name += str(count)
        return feature_name

    def readInReferenceAnnotation(self, annotation_file):
        reference_annotation = dict()
        with open(annotation_file) as f:
            for line in f:
                chromosome, source, feature, start, end, score, strand, frame, attribute = line.strip().split('\t')
                tr_id = attribute.split('transcript_id')[1].split('\"')[1]
                if len(attribute.split('gene_id')) > 1:
                    gene = attribute.split('gene_id')[1].split('\"')[1]
                else:
                    gene = None
                if feature == 'exon':
                    if tr_id not in reference_annotation:
                        reference_annotation[tr_id] = {'chromosome': chromosome, 'strand': strand, 'exons': [], 'gene': gene}
                    reference_annotation[tr_id]['exons'].append([int(start), int(end)])
        ## TAKE ADDITIONAL INFORMATION FROM EXON LISTS
        for tr_id in reference_annotation:
            t = reference_annotation[tr_id]
            t['exons'].sort(key=lambda x: x[0])
            t['tr_start'] = t['exons'][0][0]
            t['tr_end'] = t['exons'][len(t['exons'])-1][1]
            if t['strand'] == '+':
                t['tss'] = t['tr_start']
            if t['strand'] == '-':
                t['tss'] = t['tr_end']
            t['gene_length'] = t['tr_end'] - t['tr_start']
        return reference_annotation

    def createSearchWindowsFromAnnotation(self):
        ## VALUE USED TO MERGE SEARCH WINDOWS BY PROXIMITY
        join_window = 200

        temp = sorted(self.reference_annotation, key=lambda k: (
            self.reference_annotation[k]['strand'],
            self.reference_annotation[k]['chromosome'],
            self.reference_annotation[k]['tss'],
            self.reference_annotation[k]['gene'],
            k
            ))

        ## POPULATE TRANSCRIPT LIST FROM SORTED LIST;
        ## ADD SEARCH WINDOW EDGES TO ENTRIES
        transcript_list = []

        for transcript in temp:
            transcript_list.append({
                'transcript_ids': [transcript],
                'chromosome': self.reference_annotation[transcript]['chromosome'],
                'tss': self.reference_annotation[transcript]['tss'],
                'strand': self.reference_annotation[transcript]['strand'],
                'genes': [self.reference_annotation[transcript]['gene']],
                'hits': []
            })
            if self.reference_annotation[transcript]['strand'] == '+':
                transcript_list[-1]['start'] = transcript_list[-1]['tss'] - 1000
                ## MAKE SURE WINDOW END DOES NOT GO PAST TRANSCRIPT END
                end = transcript_list[-1]['tss'] + 999
                if end > self.reference_annotation[transcript]['tr_end']:
                    transcript_list[-1]['end'] = self.reference_annotation[transcript]['tr_end']
                else:
                    transcript_list[-1]['end'] = end
                #transcript_list[-1]['end'] = transcript_list[-1]['tss'] + 999
            elif self.reference_annotation[transcript]['strand'] == '-':
                ## MAKE SURE WINDOW START DOES NOT GO PAST TRANSCRIPT START
                start = transcript_list[-1]['tss'] - 999
                if start < self.reference_annotation[transcript]['tr_start']:
                    transcript_list[-1]['start'] = self.reference_annotation[transcript]['tr_end']
                else:
                    transcript_list[-1]['start'] = start
                #transcript_list[-1]['start'] = transcript_list[-1]['tss'] - 999
                transcript_list[-1]['end'] = transcript_list[-1]['tss'] + 1000

        merged_windows = []

        ## MERGE WINDOWS BASED PROXIMITY;
        ## IF WINDOWS ARE WITHIN JOIN THRESHOLD, THEY ARE MERGED;
        ## IF NOT, BUT STILL OVERLAPPING, MIDPOINT BECOMES BOUNDARY
        working_entry = transcript_list.pop(0)
        while len(transcript_list) != 0:
            next_entry = transcript_list.pop(0)
            if (working_entry['strand'] == next_entry['strand']) and (working_entry['chromosome'] == next_entry['chromosome']):
                if working_entry['tss'] + join_window >= next_entry['tss']:
                    working_entry['transcript_ids'].append(next_entry['transcript_ids'][0])
                    working_entry['genes'].append(next_entry['genes'][0])
                    working_entry['end'] = next_entry['end']
                elif working_entry['end'] >= next_entry['start']:
                    working_entry['end'] = int(math.floor((working_entry['end']+next_entry['start'])/2))
                    next_entry['start'] = working_entry['end'] + 1
                    merged_windows.append(working_entry)
                    working_entry = next_entry
                else:
                    merged_windows.append(working_entry)
                    working_entry = next_entry
            else:
                merged_windows.append(working_entry)
                working_entry = next_entry
        merged_windows.append(working_entry)

        return merged_windows

    def combineAndSortBedGraphs(self, forward_bedgraph, reverse_bedgraph):

        def readBedGraph(bedgraph_list, bedgraph_fn, strand):
            with open(bedgraph_fn) as f:
                for line in f:
                    if not ('track' in line or line == '\n'):
                        chromosome, start, end, reads = line.strip().split()
                        for i in range(int(start), int(end)):
                            bedgraph_list.append({
                                'chromosome': chromosome,
                                'start': i,
                                'end': i,
                                'reads': int(reads),
                                'strand': strand
                                })

        combined_list = []
        readBedGraph(combined_list, forward_bedgraph, '+')
        readBedGraph(combined_list, reverse_bedgraph, '-')
        return self.sortList(combined_list, 'sort_by_strand')

    ## CONSIDERS TAB-DELIMITED CHROM_SIZES FILE (UCSC)
    def findGenomeSize(self, chrom_sizes):
        genome_size = 0
        with open(chrom_sizes) as f:
            for line in f:
                genome_size += int(line.strip().split()[1])
        return genome_size

    ## FIND THRESHOLD FOR TSS CALLING, BASED ON
    ## JOTHI ET AL. (2008) NUCLEIC ACIDS RES 36: 5221-5231.
    def findReadThreshold(self, bedgraph_list, genome_size):

        def countLoci(bedgraph_list, value):
            loci = 0
            for entry in bedgraph_list:
                if entry['reads'] >= value:
                    loci += 1
            return loci

        ## IF FDR OR EXPECTED FALSE POSITIVES IS NOT SET, DEFAULT TO 1 EXPECTED FALSE POSITIVE
        if self.fdr_threshold is None and self.false_positives is None:
            self.false_positives = 1
        mappable_size = 0.8 * 2 * float(genome_size)

        read_count = 0
        for entry in bedgraph_list:
            read_count += entry['reads']
        expected_count = float(read_count)/mappable_size

        cume_probability = ((expected_count**0)/math.factorial(0))*math.exp(-expected_count)
        threshold = 1
        while True:
            probability = 1 - cume_probability
            expected_loci = probability * mappable_size
            if self.fdr_threshold:
                observed_loci = countLoci(bedgraph_list, threshold)
                fdr = float(expected_loci)/float(observed_loci)
                if fdr < fdr_threshold:
                    return threshold
            else:
                if expected_loci < self.false_positives:
                    return threshold
            cume_probability += ((expected_count**threshold)/math.factorial(threshold))*math.exp(-expected_count)
            threshold += 1

    ## FIND INTERSECTION WITH SEARCH_WINDOWS, BEDGRAPH_LIST;
    ## HITS ARE ADDED TO WINDOW_LIST, REQUIRES SORTED LIST
    def findIntersectionWithBedGraph(self, search_windows, bedgraph_list):
        search_index = 0
        bedgraph_index = 0
        while (search_index < len(search_windows)) and (bedgraph_index < len(bedgraph_list)):
            if self.isWithin(bedgraph_list[bedgraph_index], search_windows[search_index]):
                search_windows[search_index]['hits'].append([
                    bedgraph_list[bedgraph_index]['start'],
                    bedgraph_list[bedgraph_index]['reads']
                    ])
                bedgraph_index += 1
            else:
                if self.isLessThan(bedgraph_list[bedgraph_index], search_windows[search_index]):
                    bedgraph_index += 1
                else:
                    search_index += 1

    ## CALL TSS FROM ANNOTATION, REQUIRES POPULATED HITS IN WINDOWS
    def callAnnotatedTSSsFromIntersection(self, intersection, read_threshold):
        for entry in intersection:
            max_reads = float('-inf')
            max_position = None

            for hit in entry['hits']:
                if hit[1] > max_reads:
                    max_position = hit[0]
                    max_reads = hit[1]
                elif hit[1] == max_reads:
                    if entry['strand'] == '+':
                        if hit[0] < max_position:
                            max_position = hit[0]
                    elif entry['strand'] == '-':
                        if hit[0] > max_position:
                            max_position = hit[0]

            if max_reads >= read_threshold:
                self.tss_list.append({
                    'id': self.getID('obsTSS', self.annotated_tss_count),
                    'type': 'from annotation',
                    'start': max_position,
                    'end': max_position,
                    'reads': max_reads,
                    'transcript_ids': entry['transcript_ids'],
                    'genes': entry['genes'],
                    'strand': entry['strand'],
                    'chromosome': entry['chromosome']
                    })
                self.annotated_tss_count += 1

    ## CREATE WINDOWS ABOUT KNOWN TSS FOR UNANNOTATED TSSs CALLING;
    ## CONSIDERS ANNOTATED AND CALLED TSSs IN INSTANCE LISTS
    def createFilterWindowsFromAnnotationAndCalledTSSs(self):

        filter_windows = []

        if self.reference_annotation:
            for transcript in self.reference_annotation:
                filter_windows.append({
                    'strand': self.reference_annotation[transcript]['strand'],
                    'chromosome': self.reference_annotation[transcript]['chromosome'],
                    'start': self.reference_annotation[transcript]['tss'] - self.nutss_filter_size,
                    'end': self.reference_annotation[transcript]['tss'] + self.nutss_filter_size
                    })

        if self.tss_list != []:
            for tss in self.tss_list:
                filter_windows.append({
                    'strand': tss['strand'],
                    'chromosome': tss['chromosome'],
                    'start': tss['start'] - self.nutss_filter_size,
                    'end': tss['start'] + self.nutss_filter_size
                    })

        return self.sortList(filter_windows, 'sort_by_strand')

    def filterByWindowsAndThreshold(self, bedgraph_list, filter_windows, read_threshold):
        ## FILTER BY READ THRESHOLD
        filtered_bedgraph_list = []
        for entry in bedgraph_list:
            if entry['reads'] > read_threshold:
                filtered_bedgraph_list.append(entry)

        ## FILTER BY OVERLAP WITH FILTER WINDOWS
        if filter_windows != []:
            filter_index = 0
            bedgraph_index = 0
            temp = []
            while (filter_index < len(filter_windows)) and (bedgraph_index < len(filtered_bedgraph_list)):
                if self.isWithin(filtered_bedgraph_list[bedgraph_index], filter_windows[filter_index]):
                    bedgraph_index += 1
                else:
                    if self.isLessThan(filtered_bedgraph_list[bedgraph_index], filter_windows[filter_index]):
                        temp.append(filtered_bedgraph_list[bedgraph_index])
                        bedgraph_index += 1
                    else:
                        filter_index += 1
            filtered_bedgraph_list = temp
        return filtered_bedgraph_list

    ## CREATES WINDOWS FOR UNANNOTATED TSS CALLING
    def createUnannotatedSearchWindowsFromBedgraph(self, bedgraph_list):
        windows = []
        for entry in bedgraph_list:
            windows.append({
                'strand': entry['strand'],
                'chromosome': entry['chromosome'],
                'start': entry['start'] - self.nutss_search_window,
                'end': entry['end'] + self.nutss_search_window,
                'hits': []
                })

        ## MERGE OVERLAPPING WINDOWS
        merged_windows = []
        working_entry = windows.pop(0)
        while len(windows) != 0:
            next_entry = windows.pop(0)
            if (working_entry['strand'] == next_entry['strand']) and (working_entry['chromosome'] == next_entry['chromosome']):
                if working_entry['end']  >= next_entry['start']:
                    working_entry['end'] = next_entry['end']
                else:
                    merged_windows.append(working_entry)
                    working_entry = next_entry
            else:
                merged_windows.append(working_entry)
                working_entry = next_entry
        return merged_windows

    ## CALLS UNANNOTATED TSSs FROM INTERSECTION IN WINDOW HITS;
    ## PERFORMS RECURSION OVER HITS TO ENSURE THAT CALLS ARE NOT WITHIN
    ## A SEARCH WINDOW LENGTH OF EACH OTHER
    def callUnannotatedTSSsFromIntersection(self, intersection):
        for entry in intersection:
            temp = entry

            while len(temp['hits']) != 0:
                max_reads = float('-inf')
                max_position = None

                for hit in temp['hits']:
                    if hit[1] > max_reads:
                        max_position = hit[0]
                        max_reads = hit[1]
                    elif hit[1] == max_reads:
                        if temp['strand'] == '+':
                            if hit[0] < max_position:
                                max_position = hit[0]
                        elif temp['strand'] == '-':
                            if hit[0] > max_position:
                                max_position = hit[0]

                self.tss_list.append({
                    'id': self.getID('nuTSS', self.unannotated_tss_count),
                    'type': 'unannotated',
                    'start': max_position,
                    'end': max_position,
                    'reads': max_reads,
                    'strand': entry['strand'],
                    'chromosome': entry['chromosome']
                    })
                self.unannotated_tss_count += 1

                for hit in temp['hits']:
                    if abs(hit[0] - max_position) <= self.nutss_search_window:
                        temp['hits'].remove(hit)

    ## SORT CALLED TSSs AND ASSOCIATE INTO BIDIRECTIONAL PAIRS
    def associateBidirectionalTSSs(self):
        self.tss_list = self.sortList(self.tss_list, 'ignore_strand')
        for i in range(len(self.tss_list)-1):
            if self.tss_list[i]['chromosome'] == self.tss_list[i+1]['chromosome']:
                if self.tss_list[i]['strand'] == '-' and self.tss_list[i+1]['strand'] == '+':
                    if self.tss_list[i]['start'] + self.bidirectional_threshold >= self.tss_list[i+1]['start']:
                        self.tss_list[i]['partner'] = self.tss_list[i+1]['id']
                        self.tss_list[i+1]['partner'] = self.tss_list[i]['id']

    def findTSSExonIntronOverlap(self):
        exons = []
        introns = []

        if self.reference_annotation:
            for transcript in self.reference_annotation:
                for i in range(len(self.reference_annotation[transcript]['exons'])):
                    strand = self.reference_annotation[transcript]['strand']
                    chromosome = self.reference_annotation[transcript]['chromosome']
                    start = self.reference_annotation[transcript]['exons'][i][0]
                    end = self.reference_annotation[transcript]['exons'][i][1]
                    exons.append({
                        'strand': strand,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end
                        })

                for i in range(len(self.reference_annotation[transcript]['exons'])-1):
                    strand = self.reference_annotation[transcript]['strand']
                    chromosome = self.reference_annotation[transcript]['chromosome']
                    start = self.reference_annotation[transcript]['exons'][i][1]+1
                    end = self.reference_annotation[transcript]['exons'][i+1][0]-1
                    introns.append({
                        'strand': strand,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end
                        })

            exons = self.sortList(exons, 'sort_by_strand')
            introns = self.sortList(introns, 'sort_by_strand')
            self.tss_list = self.sortList(self.tss_list, 'sort_by_strand')

        def findFeatureOverlap(tss_list, feature_list, feature_key):
            if feature_list == []:
                for tss in tss_list:
                    tss[feature_key] = False
            else:
                feature_index = 0
                tss_index = 0
                while (feature_index < len(feature_list)) and (tss_index < len(tss_list)):
                    if self.isWithin(tss_list[tss_index], feature_list[feature_index]):
                        tss_list[tss_index][feature_key] = True
                        tss_index += 1
                    else:
                        if self.isLessThan(tss_list[tss_index], feature_list[feature_index]):
                            tss_list[tss_index][feature_key] = False
                            tss_index += 1
                        else:
                            feature_index += 1

        findFeatureOverlap(self.tss_list, exons, 'exon_overlap')
        findFeatureOverlap(self.tss_list, introns, 'intron_overlap')

    ## ASSOCIATE TSSs INTO CLUSTERS BY PROXIMITY;
    ## NOTE TSS CLUSTER AND NUMBER OF TSSs IN ASSOCIATED CLUSTER IN TSS ENTRY
    def associateTSSsIntoClusters(self):
        cluster_count = dict()
        self.tss_list = self.sortList(self.tss_list, 'ignore_strand')

        current_cluster = self.getID('cluster', self.tss_cluster_count)
        self.tss_cluster_count += 1
        self.tss_list[0]['cluster'] = current_cluster
        cluster_count[current_cluster] = 1

        for i in range(1, len(self.tss_list)):
            if not (self.tss_list[i-1]['chromosome'] == self.tss_list[i]['chromosome'] and self.tss_list[i-1]['start'] + self.cluster_threshold >= self.tss_list[i]['start']):
                current_cluster = self.getID('cluster', self.tss_cluster_count)
                self.tss_cluster_count += 1
            self.tss_list[i]['cluster'] = current_cluster
            if current_cluster not in cluster_count:
                cluster_count[current_cluster] = 1
            else:
                cluster_count[current_cluster] += 1

        for tss in self.tss_list:
            tss['cluster_count'] = cluster_count[tss['cluster']]

    def createDetailFile(self):
        self.findTSSExonIntronOverlap()
        self.associateTSSsIntoClusters()
        with open(self.detail_file, 'w') as OUTPUT:
            OUTPUT.write('TSS ID\tType\tTranscripts\tGenes\tStrand\tChromosome\tPosition\tReads\tBidirectional?\tBidirectional partner\tGene overlap?\tExon overlap?\tIntron overlap?\tTSS cluster\tTSSs in associated cluster\n')
            for tss in self.tss_list:
                OUTPUT.write(tss['id'])
                OUTPUT.write('\t' + tss['type'])
                for entry in ['transcript_ids', 'genes']:
                    if entry in tss:
                        OUTPUT.write('\t' + ','.join(tss[entry]))
                    else:
                        OUTPUT.write('\tNA')
                for entry in ['strand', 'chromosome', 'start', 'reads']:
                    OUTPUT.write('\t' + str(tss[entry]))
                if 'partner' in tss:
                    OUTPUT.write('\tTrue\t' + tss['partner'])
                else:
                    OUTPUT.write('\tFalse\tNA')
                OUTPUT.write('\t' + str(tss['exon_overlap'] or tss['intron_overlap']))
                for entry in ['exon_overlap', 'intron_overlap', 'cluster', 'cluster_count']:
                    OUTPUT.write('\t' + str(tss[entry]))
                OUTPUT.write('\n')

    def writeBedFile(self, tss_list, output_bed):
        with open(output_bed, 'w') as OUTPUT:
            for tss in tss_list:
                OUTPUT.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    tss['chromosome'],
                    str(tss['start'] - 1),
                    str(tss['start']),
                    tss['id'],
                    '0',
                    tss['strand']
                    ))

    def callTSSsFromAnnotation(self, bedgraph_list, read_threshold):
        ref_search_windows = self.createSearchWindowsFromAnnotation()
        self.findIntersectionWithBedGraph(ref_search_windows, bedgraph_list)
        self.callAnnotatedTSSsFromIntersection(ref_search_windows, read_threshold)

    def callUnannotatedTSSs(self, bedgraph_list, read_threshold):
        filter_windows = self.createFilterWindowsFromAnnotationAndCalledTSSs()
        filtered_bedgraph = self.filterByWindowsAndThreshold(bedgraph_list, filter_windows, read_threshold)
        unannotated_search_windows = self.createUnannotatedSearchWindowsFromBedgraph(filtered_bedgraph)
        self.findIntersectionWithBedGraph(unannotated_search_windows, filtered_bedgraph)
        self.callUnannotatedTSSsFromIntersection(unannotated_search_windows)

    def execute(self):
        bedgraph_list = self.combineAndSortBedGraphs(self.forward_bedgraph, self.reverse_bedgraph)
        genome_size = self.findGenomeSize(self.chrom_sizes)
        read_threshold = self.findReadThreshold(bedgraph_list, genome_size)

        if self.annotation_file:
            self.reference_annotation = self.readInReferenceAnnotation(self.annotation_file)
            self.callTSSsFromAnnotation(bedgraph_list, read_threshold)
        self.callUnannotatedTSSs(bedgraph_list, read_threshold)
        self.associateBidirectionalTSSs()
        if self.detail_file:
            self.createDetailFile()
        self.writeBedFile(self.tss_list, self.output_bed)

parser = argparse.ArgumentParser()
parser.add_argument('--fdr', default=None, type=float, help='set read threshold by FDR (FLOAT)')
parser.add_argument('--false_positives', default=None, type=int, help='set read threshold by false positive count (INTEGER)')
parser.add_argument('--nutss_filter_size', default=750, type=int, help='set nuTSS filter size; any read within INTEGER of obsTSS/annoTSS is filtered prior to nuTSS calling')
parser.add_argument('--nutss_search_window', default=250, type=int, help='set nuTSS search window size to INTEGER')
parser.add_argument('--bidirectional_threshold', default=1000, type=int, help='INTEGER threshold to associate bidirectional TSSs')
parser.add_argument('--detail_file', default=None, type=str, help='create a tab-delimited TXT file with details about TSS calls')
parser.add_argument('--cluster_threshold', default=1000, type=int, help='INTEGER threshold to associate TSSs into clusters')
parser.add_argument('--annotation_file', '-a', type=str, help='annotation in GTF format')
parser.add_argument('forward_bedgraph', type=str, help='forward strand Start-seq bedgraph file')
parser.add_argument('reverse_bedgraph', type=str, help='reverse strand Start-seq bedgraph file')
parser.add_argument('chrom_sizes', type=str, help='standard tab-delimited chromosome sizes file')
parser.add_argument('output_bed', type=str, help='output TSS BED file')
args = parser.parse_args()

TSSCalling(**vars(args))
