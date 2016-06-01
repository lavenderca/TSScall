#!/usr/bin/env python

# CREATED BY CHRISTOPHER LAVENDER
# BASED ON WORK BY ADAM BURKHOLDER
# INTEGRATIVE BIOINFORMATICS, NIEHS
# WORKING OBJECT ORIENTED VERSION

import os
import math
import argparse
import sys
from operator import itemgetter


# STRAND_STATUS IS USED TO DETERMINE IF STRAND IS USED IN SORT
def sortList(input_list, strand_status):
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


# ENTRY 1 IS LESS THAN ENTRY 2?
def isLessThan(entry_1, entry_2):
    for val in ['strand', 'chromosome', 'start']:
        if entry_1[val] < entry_2[val]:
            return True
        elif entry_1[val] > entry_2[val]:
            return False
    return False


# ENTRY 1 IS WITHIN ENTRY 2?
def isWithin(entry_1, entry_2):
    if entry_1['strand'] == entry_2['strand'] and\
            entry_1['chromosome'] == entry_2['chromosome']:
        if entry_1['start'] >= entry_2['start'] and\
                entry_1['end'] <= entry_2['end']:
            return True
    return False


def getID(base_name, count):
    max_entries = 999999
    feature_name = base_name + '_'
    for i in range(len(str(count)), len(str(max_entries))):
        feature_name += '0'
    feature_name += str(count)
    return feature_name


def readInReferenceAnnotation(annotation_file):
    reference_annotation = dict()
    all_gtf_keys = []
    with open(annotation_file) as f:
        for line in f:
            chromosome, source, feature, start, end, score, strand, frame,\
                attributes = line.strip().split('\t')

            keys = []
            values = []
            gtf_fields = dict()
            for entry in attributes.split(';')[:-1]:
                keys.append(entry.split('\"')[0].strip())
                values.append(entry.split('\"')[1].strip())
            for key, value in zip(keys, values):
                gtf_fields[key] = [value]
            for key in keys:
                if key not in all_gtf_keys:
                    all_gtf_keys.append(key)

            tr_id = gtf_fields.pop('transcript_id')[0]
            gene_id = gtf_fields.pop('gene_id')[0]
            for val in ('transcript_id', 'gene_id'):
                all_gtf_keys.remove(val)

            if feature == 'exon':
                ref_id = (tr_id, chromosome)
                if ref_id not in reference_annotation:
                    reference_annotation[ref_id] = {
                        'chromosome': chromosome,
                        'strand': strand,
                        'exons': [],
                        'gene_id': gene_id,
                        'gtf_fields': gtf_fields,
                        }
                reference_annotation[ref_id]['exons'].append(
                    [int(start), int(end)]
                    )
    for ref_id in reference_annotation:
        t = reference_annotation[ref_id]
        # TAKE ADDITIONAL INFORMATION FROM EXON LISTS
        t['exons'].sort(key=lambda x: x[0])
        t['tr_start'] = t['exons'][0][0]
        t['tr_end'] = t['exons'][len(t['exons'])-1][1]
        if t['strand'] == '+':
            t['tss'] = t['tr_start']
        if t['strand'] == '-':
            t['tss'] = t['tr_end']
        t['gene_length'] = t['tr_end'] - t['tr_start']
        # POPULATE MISSING GTF FIELD ENTRIES
        for key in all_gtf_keys:
            if key not in t['gtf_fields']:
                t['gtf_fields'][key] = [None]
    return reference_annotation, all_gtf_keys


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
        self.cluster_bed = kwargs['cluster_bed']
        self.call_method = kwargs['call_method']
        self.annotation_join_distance = kwargs['annotation_join_distance']
        self.annotation_search_window = kwargs['annotation_search_window']
        self.bin_winner_size = kwargs['bin_winner_size']

        self.set_read_threshold = kwargs['set_read_threshold']
        try:
            int(self.set_read_threshold)
        except:
            pass
        else:
            self.set_read_threshold = int(self.set_read_threshold)

        # EVALUATE THRESHOLD METHOD ARGUMENTS; IF NONE, SET FDR_THRESHOLD
        # AT 0.001
        implied_threshold_methods = 0
        for val in [
                self.fdr_threshold,
                self.false_positives,
                self.set_read_threshold]:
            implied_threshold_methods += int(bool(val))
        if implied_threshold_methods == 1:
            pass
        elif implied_threshold_methods > 1:
            raise ValueError('More than 1 read threshold method implied!!')
        elif implied_threshold_methods == 0:
            self.fdr_threshold = 0.001

        self.tss_list = []
        self.reference_annotation = None
        self.gtf_attribute_fields = []
        self.annotated_tss_count = 0
        self.unannotated_tss_count = 0
        self.tss_cluster_count = 0
        self.unobserved_ref_count = 0

        self.execute()

    def createSearchWindowsFromAnnotation(self):
        # VALUE USED TO MERGE SEARCH WINDOWS BY PROXIMITY
        join_window = self.annotation_join_distance
        window_size = self.annotation_search_window

        current_entry = sorted(self.reference_annotation, key=lambda k: (
            self.reference_annotation[k]['strand'],
            self.reference_annotation[k]['chromosome'],
            self.reference_annotation[k]['tss'],
            # self.reference_annotation[k]['gene'],
            k,
            ))

        # POPULATE TRANSCRIPT LIST FROM SORTED LIST;
        # ADD SEARCH WINDOW EDGES TO ENTRIES
        transcript_list = []

        for ref in current_entry:
            transcript_list.append({
                'transcript_id': [ref[0]],
                'chromosome':
                    self.reference_annotation[ref]['chromosome'],
                'tss': [self.reference_annotation[ref]['tss']],
                'strand': self.reference_annotation[ref]['strand'],
                'gene_id': [self.reference_annotation[ref]['gene_id']],
                'hits': [],
                'gtf_fields': self.reference_annotation[ref]['gtf_fields'],
            })
            if self.reference_annotation[ref]['strand'] == '+':
                transcript_list[-1]['start'] = \
                    transcript_list[-1]['tss'][0] - window_size
                # MAKE SURE WINDOW END DOES NOT GO PAST TRANSCRIPT END
                end = transcript_list[-1]['tss'][0] + window_size
                if end > self.reference_annotation[ref]['tr_end']:
                    transcript_list[-1]['end'] = \
                        self.reference_annotation[ref]['tr_end']
                else:
                    transcript_list[-1]['end'] = end
            elif self.reference_annotation[ref]['strand'] == '-':
                # MAKE SURE WINDOW START DOES NOT GO PAST TRANSCRIPT START
                start = transcript_list[-1]['tss'][0] - window_size
                if start < self.reference_annotation[ref]['tr_start']:
                    transcript_list[-1]['start'] = \
                        self.reference_annotation[ref]['tr_end']
                else:
                    transcript_list[-1]['start'] = start
                transcript_list[-1]['end'] = \
                    transcript_list[-1]['tss'][0] + window_size

        merged_windows = []

        # MERGE WINDOWS BASED PROXIMITY;
        # IF WINDOWS ARE WITHIN JOIN THRESHOLD, THEY ARE MERGED;
        # IF NOT, BUT STILL OVERLAPPING, MIDPOINT BECOMES BOUNDARY
        working_entry = transcript_list.pop(0)
        while len(transcript_list) != 0:
            next_entry = transcript_list.pop(0)
            if (working_entry['strand'] == next_entry['strand']) and \
                    (working_entry['chromosome'] == next_entry['chromosome']):
                if working_entry['tss'][-1] + join_window >= \
                        next_entry['tss'][0]:
                    working_entry['transcript_id'].append(
                        next_entry['transcript_id'][0]
                        )
                    working_entry['gene_id'].append(
                        next_entry['gene_id'][0]
                        )
                    for key in working_entry['gtf_fields']:
                        working_entry['gtf_fields'][key].append(
                            next_entry['gtf_fields'][key][0]
                            )
                    # working_entry['genes'].append(next_entry['genes'][0])
                    working_entry['end'] = next_entry['end']
                    working_entry['tss'].append(next_entry['tss'][0])
                elif working_entry['end'] >= next_entry['start']:
                    working_entry['end'] = int(math.floor(
                        (working_entry['end']+next_entry['start'])/2
                        ))
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
                        for i in range(int(start)+1, int(end)+1):
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
        return sortList(combined_list, 'sort_by_strand')

    # CONSIDERS TAB-DELIMITED CHROM_SIZES FILE (UCSC)
    def findGenomeSize(self, chrom_sizes):
        genome_size = 0
        with open(chrom_sizes) as f:
            for line in f:
                genome_size += int(line.strip().split()[1])
        return genome_size

    # FIND THRESHOLD FOR TSS CALLING, BASED ON
    # JOTHI ET AL. (2008) NUCLEIC ACIDS RES 36: 5221-5231.
    def findReadThreshold(self, bedgraph_list, genome_size):

        def countLoci(bedgraph_list, value):
            loci = 0
            for entry in bedgraph_list:
                if entry['reads'] >= value:
                    loci += 1
            return loci

        if self.fdr_threshold or self.false_positives:
            self.false_positives = 1
            mappable_size = 0.8 * 2 * float(genome_size)
            read_count = 0
            for entry in bedgraph_list:
                read_count += entry['reads']
            expected_count = float(read_count)/mappable_size

            cume_probability = ((expected_count**0)/math.factorial(0)) * \
                math.exp(-expected_count)
            threshold = 1
            while True:
                probability = 1 - cume_probability
                expected_loci = probability * mappable_size
                if self.fdr_threshold:
                    observed_loci = countLoci(bedgraph_list, threshold)
                    fdr = float(expected_loci)/float(observed_loci)
                    if fdr < self.fdr_threshold:
                        return threshold
                else:
                    if expected_loci < self.false_positives:
                        return threshold
                cume_probability += \
                    ((expected_count**threshold)/math.factorial(threshold)) * \
                    math.exp(-expected_count)
                threshold += 1
        else:
            return self.set_read_threshold

    # FIND INTERSECTION WITH SEARCH_WINDOWS, BEDGRAPH_LIST;
    # HITS ARE ADDED TO WINDOW_LIST, REQUIRES SORTED LIST
    def findIntersectionWithBedGraph(self, search_windows, bedgraph_list):
        search_index = 0
        bedgraph_index = 0
        while (search_index < len(search_windows)) and \
                (bedgraph_index < len(bedgraph_list)):
            if isWithin(bedgraph_list[bedgraph_index],
                        search_windows[search_index]):
                search_windows[search_index]['hits'].append([
                    bedgraph_list[bedgraph_index]['start'],
                    bedgraph_list[bedgraph_index]['reads']
                    ])
                bedgraph_index += 1
            else:
                if isLessThan(bedgraph_list[bedgraph_index],
                              search_windows[search_index]):
                    bedgraph_index += 1
                else:
                    search_index += 1

    # CREATE WINDOWS ABOUT KNOWN TSS FOR UNANNOTATED TSSs CALLING;
    # CONSIDERS ANNOTATED AND CALLED TSSs IN INSTANCE LISTS
    def createFilterWindowsFromAnnotationAndCalledTSSs(self):

        filter_windows = []

        if self.reference_annotation:
            for transcript in self.reference_annotation:
                filter_windows.append({
                    'strand': self.reference_annotation[transcript]['strand'],
                    'chromosome':
                        self.reference_annotation[transcript]['chromosome'],
                    'start': self.reference_annotation[transcript]['tss'] -
                        self.nutss_filter_size,
                    'end': self.reference_annotation[transcript]['tss'] +
                        self.nutss_filter_size
                    })

        if self.tss_list != []:
            for tss in self.tss_list:
                filter_windows.append({
                    'strand': tss['strand'],
                    'chromosome': tss['chromosome'],
                    'start': tss['start'] - self.nutss_filter_size,
                    'end': tss['start'] + self.nutss_filter_size
                    })

        return sortList(filter_windows, 'sort_by_strand')

    def filterBedGraphListByWindows(self, bedgraph_list, filter_windows):
        # FILTER BY OVERLAP WITH FILTER WINDOWS
        if filter_windows != []:
            filter_index = 0
            bedgraph_index = 0
            working_list = []
            while (filter_index < len(filter_windows)) and \
                    (bedgraph_index < len(bedgraph_list)):
                if isWithin(bedgraph_list[bedgraph_index],
                            filter_windows[filter_index]):
                    bedgraph_index += 1
                else:
                    if isLessThan(bedgraph_list[bedgraph_index],
                                  filter_windows[filter_index]):
                        working_list.append(bedgraph_list[bedgraph_index])
                        bedgraph_index += 1
                    else:
                        filter_index += 1
            bedgraph_list = working_list
        return bedgraph_list

    # CREATES WINDOWS FOR UNANNOTATED TSS CALLING
    def createUnannotatedSearchWindowsFromBedgraph(self,
                                                   bedgraph_list,
                                                   read_threshold):
        windows = []
        for entry in bedgraph_list:
            if entry['reads'] > read_threshold:
                windows.append({
                    'strand': entry['strand'],
                    'chromosome': entry['chromosome'],
                    'start': entry['start'] - self.nutss_search_window,
                    'end': entry['end'] + self.nutss_search_window,
                    'hits': []
                    })

        # MERGE OVERLAPPING WINDOWS
        merged_windows = []
        working_entry = windows.pop(0)
        while len(windows) != 0:
            next_entry = windows.pop(0)
            if (working_entry['strand'] == next_entry['strand']) and\
                    (working_entry['chromosome'] == next_entry['chromosome']):
                if working_entry['end'] >= next_entry['start']:
                    working_entry['end'] = next_entry['end']
                else:
                    merged_windows.append(working_entry)
                    working_entry = next_entry
            else:
                merged_windows.append(working_entry)
                working_entry = next_entry
        return merged_windows

    # SORT CALLED TSSs AND ASSOCIATE INTO BIDIRECTIONAL PAIRS
    def associateBidirectionalTSSs(self):
        self.tss_list = sortList(self.tss_list, 'ignore_strand')
        for i in range(len(self.tss_list)-1):
            if self.tss_list[i]['chromosome'] == \
                    self.tss_list[i+1]['chromosome']:
                if self.tss_list[i]['strand'] == '-' and \
                        self.tss_list[i+1]['strand'] == '+':
                    if self.tss_list[i]['start'] + \
                            self.bidirectional_threshold >= \
                            self.tss_list[i+1]['start']:
                        distance = abs(self.tss_list[i]['start'] -
                                       self.tss_list[i+1]['start'])
                        self.tss_list[i]['divergent partner'] = \
                            self.tss_list[i+1]['id']
                        self.tss_list[i+1]['divergent partner'] = \
                            self.tss_list[i]['id']
                        self.tss_list[i]['divergent distance'] = distance
                        self.tss_list[i+1]['divergent distance'] = distance
                if self.tss_list[i]['strand'] == '+' and \
                        self.tss_list[i+1]['strand'] == '-':
                    if self.tss_list[i]['start'] + \
                            self.bidirectional_threshold >= \
                            self.tss_list[i+1]['start']:
                        distance = abs(self.tss_list[i]['start'] -
                                       self.tss_list[i+1]['start'])
                        self.tss_list[i]['convergent partner'] = \
                            self.tss_list[i+1]['id']
                        self.tss_list[i+1]['convergent partner'] = \
                            self.tss_list[i]['id']
                        self.tss_list[i]['convergent distance'] = distance
                        self.tss_list[i+1]['convergent distance'] = distance

    def findTSSExonIntronOverlap(self):
        exons = []
        introns = []

        if self.reference_annotation:
            for transcript in self.reference_annotation:
                for i in range(len(
                        self.reference_annotation[transcript]['exons'])):
                    strand = self.reference_annotation[transcript]['strand']
                    chromosome =\
                        self.reference_annotation[transcript]['chromosome']
                    start =\
                        self.reference_annotation[transcript]['exons'][i][0]
                    end = self.reference_annotation[transcript]['exons'][i][1]
                    exons.append({
                        'strand': strand,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end
                        })

                for i in range(
                        len(self.reference_annotation[transcript]['exons'])-1):
                    strand = self.reference_annotation[transcript]['strand']
                    chromosome =\
                        self.reference_annotation[transcript]['chromosome']
                    start = \
                        self.reference_annotation[transcript]['exons'][i][1]+1
                    end = \
                        self.reference_annotation[transcript]['exons'][i+1][0]\
                        - 1
                    introns.append({
                        'strand': strand,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end
                        })

            exons = sortList(exons, 'sort_by_strand')
            introns = sortList(introns, 'sort_by_strand')
            self.tss_list = sortList(self.tss_list, 'sort_by_strand')

        def findFeatureOverlap(tss_list, feature_list, feature_key):
            if feature_list == []:
                for tss in tss_list:
                    tss[feature_key] = False
            else:
                feature_index = 0
                tss_index = 0
                while (feature_index < len(feature_list)) and\
                        (tss_index < len(tss_list)):
                    if isWithin(tss_list[tss_index],
                                feature_list[feature_index]):
                        tss_list[tss_index][feature_key] = True
                        tss_index += 1
                    else:
                        if isLessThan(tss_list[tss_index],
                                      feature_list[feature_index]):
                            tss_list[tss_index][feature_key] = False
                            tss_index += 1
                        else:
                            feature_index += 1

        findFeatureOverlap(self.tss_list, exons, 'exon_overlap')
        findFeatureOverlap(self.tss_list, introns, 'intron_overlap')

    # ASSOCIATE TSSs INTO CLUSTERS BY PROXIMITY;
    # ADD TSS CLUSTER AND NUMBER OF TSSs IN ASSOCIATED CLUSTER IN TSS ENTRY
    def associateTSSsIntoClusters(self):
        cluster_count = dict()
        self.tss_list = sortList(self.tss_list, 'ignore_strand')

        current_cluster = getID('cluster', self.tss_cluster_count)
        self.tss_cluster_count += 1
        self.tss_list[0]['cluster'] = current_cluster
        cluster_count[current_cluster] = 1

        for i in range(1, len(self.tss_list)):
            if not (self.tss_list[i-1]['chromosome'] ==
                    self.tss_list[i]['chromosome'] and
                    self.tss_list[i-1]['start'] + self.cluster_threshold >=
                    self.tss_list[i]['start']):
                current_cluster = getID('cluster', self.tss_cluster_count)
                self.tss_cluster_count += 1
            self.tss_list[i]['cluster'] = current_cluster
            if current_cluster not in cluster_count:
                cluster_count[current_cluster] = 1
            else:
                cluster_count[current_cluster] += 1

        for tss in self.tss_list:
            tss['cluster_count'] = cluster_count[tss['cluster']]

    def createDetailFile(self):
        def checkHits(window):
            for hit in window['hits']:
                if hit[1] >= self.read_threshold:
                    return True
            return False

        def writeUnobservedEntry(OUTPUT, tss, tr_ids, gene_ids, window):
            tss_id = getID('annoTSS', self.unobserved_ref_count)
            self.unobserved_ref_count += 1

            transcripts = tr_ids[0]
            genes = gene_ids[0]
            for i in range(1, len(tr_ids)):
                transcripts += ';' + tr_ids[i]
                genes += ';' + gene_ids[i]

            reads = 0
            for hit in window['hits']:
                if int(tss) == int(hit[0]):
                    reads = hit[1]

            OUTPUT.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\
                         \t{}\t{}\t{}'
                         .format(
                            tss_id,
                            'unobserved reference TSS',
                            transcripts,
                            genes,
                            window['strand'],
                            window['chromosome'],
                            str(tss),
                            str(reads),
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                         ))
            for key in self.gtf_attribute_fields:
                # OUTPUT.write('\t' + ';'.join(window['gtf_fields'][key]))
                OUTPUT.write('\t' + ';'.join(['None' if v is None else v for
                                              v in window['gtf_fields'][key]]))
            OUTPUT.write('\n')

        # self.findTSSExonIntronOverlap()
        # self.associateTSSsIntoClusters()
        # Remove GTF fields 'exon_number' and 'exon_id' if present
        skip_fields = ['exon_number', 'exon_id']
        for entry in skip_fields:
            if entry in self.gtf_attribute_fields:
                self.gtf_attribute_fields.remove(entry)

        with open(self.detail_file, 'w') as OUTPUT:
            OUTPUT.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\
                         \t{}\t{}\t{}'
                         .format(
                            'TSS ID',
                            'Type',
                            'Transcripts',
                            'Gene ID',
                            'Strand',
                            'Chromsome',
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
                         ))
            for field in self.gtf_attribute_fields:
                OUTPUT.write('\t' + field)
            OUTPUT.write('\n')

            for tss in self.tss_list:
                OUTPUT.write(tss['id'])
                OUTPUT.write('\t' + tss['type'])
                for key in ('transcript_id', 'gene_id'):
                    if key in tss:
                        OUTPUT.write('\t' + ';'.join(tss[key]))
                    else:
                        OUTPUT.write('\tNA')
                for entry in ['strand', 'chromosome', 'start', 'reads']:
                    OUTPUT.write('\t' + str(tss[entry]))
                if 'divergent partner' in tss:
                    OUTPUT.write('\t{}\t{}\t{}'.format(
                        'True',
                        tss['divergent partner'],
                        str(tss['divergent distance']),
                    ))
                else:
                    OUTPUT.write('\tFalse\tNA\tNA')
                if 'convergent partner' in tss:
                    OUTPUT.write('\t{}\t{}\t{}'.format(
                        'True',
                        tss['convergent partner'],
                        str(tss['convergent distance']),
                    ))
                else:
                    OUTPUT.write('\tFalse\tNA\tNA')
                # OUTPUT.write('\t' + str(
                #     tss['exon_overlap'] or tss['intron_overlap']))
                for entry in [
                        'cluster',
                        'cluster_count']:
                    OUTPUT.write('\t' + str(tss[entry]))
                if 'gtf_fields' in tss:
                    for key in self.gtf_attribute_fields:
                        # OUTPUT.write('\t' + ';'.join(tss['gtf_fields'][key]))
                        OUTPUT.write('\t' + ';'.join(
                            ['None' if v is None else
                             v for v in tss['gtf_fields'][key]]
                            ))
                else:
                    for key in self.gtf_attribute_fields:
                        OUTPUT.write('\tNA')
                OUTPUT.write('\n')

            if self.annotation_file:
                for window in self.ref_search_windows:
                    if not checkHits(window):
                        window_tss = []
                        for tr_id, gene_id, tss in zip(window['transcript_id'],
                                                       window['gene_id'],
                                                       window['tss']):
                            window_tss.append({
                                'transcript_id': tr_id,
                                'gene_id': gene_id,
                                'tss': int(tss),
                            })
                        window_tss.sort(key=itemgetter('tss'))

                        current_tss = window_tss[0]['tss']
                        current_tr_ids = [window_tss[0]['transcript_id']]
                        current_genes = [window_tss[0]['gene_id']]
                        window_index = 1
                        while window_index < len(window_tss):
                            if current_tss == window_tss[window_index]['tss']:
                                current_tr_ids.append(
                                    window_tss[window_index]['transcript_id'])
                                current_genes.append(
                                    window_tss[window_index]['gene_id'])
                            else:
                                writeUnobservedEntry(OUTPUT, current_tss,
                                                     current_tr_ids,
                                                     current_genes,
                                                     window)
                                current_tss = window_tss[window_index]['tss']
                                current_tr_ids = \
                                    [window_tss[window_index]['transcript_id']]
                                current_genes = [window_tss[0]['gene_id']]
                            window_index += 1
                        writeUnobservedEntry(OUTPUT, current_tss,
                                             current_tr_ids, current_genes,
                                             window)

    def writeClusterBed(self, tss_list, cluster_bed):
        clusters = dict()
        with open(cluster_bed, 'w') as OUTPUT:
            for tss in tss_list:
                if tss['cluster'] in clusters:
                    clusters[tss['cluster']]['tss'].append(tss['start'])
                else:
                    clusters[tss['cluster']] = {
                        'chromosome': tss['chromosome'],
                        'tss': [tss['start']],
                    }
            for cluster in sorted(clusters):
                tss = sorted(clusters[cluster]['tss'])
                OUTPUT.write('{}\t{}\t{}\t{}\n'.format(
                    clusters[cluster]['chromosome'],
                    str(tss[0] - 1),
                    str(tss[-1]),
                    cluster,
                ))

    def writeBedFile(self, tss_list, output_bed):
        with open(output_bed, 'w') as OUTPUT:
            for tss in tss_list:
                OUTPUT.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    tss['chromosome'],
                    str(tss['start'] - 1),
                    str(tss['start']),
                    tss['id'],
                    '0',
                    tss['strand']
                    ))

    # FROM HITS IN SEARCH WINDOWS, CALL TSSs
    # COUNT IS RETURNED IN ORDER TO UPDATE INSTANCE VARIABLES
    def callTSSsFromIntersection(self, intersection, read_threshold, base_name,
                                 count, tss_type, nearest_allowed):
        def callTSS(hits, strand):
            if self.call_method == 'global':
                max_reads = float('-inf')
                max_position = None
                for hit in hits:
                    if hit[1] > max_reads:
                        max_position = hit[0]
                        max_reads = hit[1]
                    elif hit[1] == max_reads:
                        if strand == '+':
                            if hit[0] < max_position:
                                max_position = hit[0]
                        elif strand == '-':
                            if hit[0] > max_position:
                                max_position = hit[0]
                return max_position, max_reads
            if self.call_method == 'bin_winner':
                bin_size = self.bin_winner_size
                bins = []
                # MAKE BINS
                hits.sort(key=itemgetter(0))
                for i in range(len(hits)):
                    bins.append({
                        'total_reads': 0,
                        'bin_hits': []
                        })
                    for j in range(i, len(hits)):
                        if abs(hits[i][0] - hits[j][0]) <= bin_size:
                            bins[-1]['total_reads'] += hits[j][1]
                            bins[-1]['bin_hits'].append(hits[j])
                # SELECT BIN WITH HIGHEST TOTAL READS
                # BECAUSE SORTED, WILL TAKE UPSTREAM BIN IN TIES
                max_bin_reads = float('-inf')
                max_bin_index = None
                for i, entry in enumerate(bins):
                    if entry['total_reads'] > max_bin_reads:
                        max_bin_index = i
                        max_bin_reads = entry['total_reads']
                # GET LOCAL WINNER
                # BECAUSE SORTED, WILL TAKE UPSTREAM TSS IN TIES
                max_reads = float('-inf')
                max_position = None
                for hit in bins[max_bin_index]['bin_hits']:
                    if hit[1] > max_reads:
                        max_position = hit[0]
                        max_reads = hit[1]
                return max_position, max_reads

        # ITERATE THROUGH WINDOWS IN INTERSECTION
        for entry in intersection:
            entry_hits = entry['hits']
            # LOOP WHILE 'HITS' IS POPULATED
            while len(entry_hits) != 0:
                # CALL A TSS
                tss_position, tss_reads = callTSS(entry_hits, entry['strand'])
                if tss_reads >= read_threshold:
                    self.tss_list.append({
                        'id': getID(base_name, count),
                        'type': tss_type,
                        'start': tss_position,
                        'end': tss_position,
                        'reads': tss_reads,
                        })
                    # IF VAL IN ENTRY, ADD TO DICT IN TSS LIST
                    for val in ['transcript_id', 'gene_id', 'strand',
                                'chromosome', 'gtf_fields']:
                        if val in entry:
                            self.tss_list[-1][val] = entry[val]
                    count += 1
                # GO THROUGH HITS, KEEP THOSE WITHIN NEAREST_ALLOWED
                temp = []
                for hit in entry_hits:
                    if abs(hit[0] - tss_position) > nearest_allowed:
                        temp.append(hit)
                entry_hits = temp
        return count

    def callTSSsFromAnnotation(self, bedgraph_list, read_threshold):
        self.ref_search_windows = self.createSearchWindowsFromAnnotation()
        self.findIntersectionWithBedGraph(self.ref_search_windows,
                                          bedgraph_list)
        self.annotated_tss_count = self.callTSSsFromIntersection(
            self.ref_search_windows,
            read_threshold,
            'obsTSS',
            self.annotated_tss_count,
            'called from reference window',
            float('inf')
            )

    def callUnannotatedTSSs(self, bedgraph_list, read_threshold):
        filter_windows = self.createFilterWindowsFromAnnotationAndCalledTSSs()
        filtered_bedgraph = self.filterBedGraphListByWindows(bedgraph_list,
                                                             filter_windows)
        unannotated_search_windows =\
            self.createUnannotatedSearchWindowsFromBedgraph(filtered_bedgraph,
                                                            read_threshold)
        self.findIntersectionWithBedGraph(unannotated_search_windows,
                                          filtered_bedgraph)
        self.unannotated_tss_count = self.callTSSsFromIntersection(
            unannotated_search_windows,
            read_threshold,
            'nuTSS',
            self.unannotated_tss_count,
            'unannotated',
            self.nutss_search_window
            )

    def execute(self):
        sys.stdout.write('Reading in bedGraph files...\n')
        bedgraph_list = self.combineAndSortBedGraphs(self.forward_bedgraph,
                                                     self.reverse_bedgraph)
        genome_size = self.findGenomeSize(self.chrom_sizes)
        sys.stdout.write('Calculating read threshold...\n')
        self.read_threshold = \
            self.findReadThreshold(bedgraph_list, genome_size)
        sys.stdout.write('Read threshold set to {}\n'.format(
            str(self.read_threshold)))

        if self.annotation_file:
            sys.stdout.write('Reading in annotation file...\n')
            self.reference_annotation, self.gtf_attribute_fields =\
                readInReferenceAnnotation(self.annotation_file)
            sys.stdout.write('Calling TSSs from annotation...\n')
            self.callTSSsFromAnnotation(bedgraph_list, self.read_threshold)
            sys.stdout.write('{} TSSs called from annotation\n'.format(
                str(self.annotated_tss_count)))
        sys.stdout.write('Calling unannotated TSSs...\n')
        self.callUnannotatedTSSs(bedgraph_list, self.read_threshold)
        sys.stdout.write('{} unannotated TSSs called\n'.format(
            str(self.unannotated_tss_count)))
        sys.stdout.write('Associating bidirectional TSSs...\n')
        self.associateBidirectionalTSSs()
        self.associateTSSsIntoClusters()
        if self.detail_file:
            sys.stdout.write('Creating detail file...\n')
            self.createDetailFile()
        if self.cluster_bed:
            sys.stdout.write('Creating cluster bed...\n')
            self.writeClusterBed(self.tss_list, self.cluster_bed)
        sys.stdout.write('Creating output bed...\n')
        self.writeBedFile(self.tss_list, self.output_bed)
        sys.stdout.write('TSS calling complete\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fdr', default=None, type=float,
                        help='set read threshold by FDR (FLOAT) (Default \
                        method: less than 0.001)')
    parser.add_argument('--false_positives', default=None, type=int,
                        help='set read threshold by false positive count')
    parser.add_argument('--nutss_filter_size', default=750, type=int,
                        help='set nuTSS filter size; any read within INTEGER \
                        of obsTSS/annoTSS is filtered prior to nuTSS calling \
                        (Default: 750)')
    parser.add_argument('--nutss_search_window', default=250, type=int,
                        help='set nuTSS search window size to INTEGER \
                        (Default: 250)')
    parser.add_argument('--bidirectional_threshold', default=1000, type=int,
                        help='INTEGER threshold to associate bidirectional \
                        TSSs (Default: 1000)')
    parser.add_argument('--detail_file', default=None, type=str,
                        help='create a tab-delimited TXT file with details \
                        about TSS calls')
    parser.add_argument('--cluster_threshold', default=1000, type=int,
                        help='INTEGER threshold to associate TSSs into \
                        clusters (Default: 1000)')
    parser.add_argument('--annotation_file', '-a', type=str,
                        help='annotation in GTF format')
    parser.add_argument('--call_method', type=str, default='bin_winner',
                        choices=['global', 'bin_winner'],
                        help='TSS calling method to use (Default: bin_winner)')
    parser.add_argument('--annotation_join_distance', type=int, default=200,
                        help='set INTEGER distace threshold for joining search \
                        windows from annotation (Default: 200)')
    parser.add_argument('--annotation_search_window', type=int, default=1000,
                        help='set annotation search window size to INTEGER \
                        (Default: 1000)')
    parser.add_argument('--set_read_threshold', type=float, default=None,
                        help='set read threshold for TSS calling to FLOAT; do \
                        not determine threshold from data')
    parser.add_argument('--bin_winner_size', type=int, default=200,
                        help='set bin size for call method bin_winner \
                        (Default: 200)')
    parser.add_argument('--cluster_bed', type=str, default=None,
                        help='write clusters to output bed file')
    parser.add_argument('forward_bedgraph', type=str,
                        help='forward strand Start-seq bedgraph file')
    parser.add_argument('reverse_bedgraph', type=str,
                        help='reverse strand Start-seq bedgraph file')
    parser.add_argument('chrom_sizes', type=str,
                        help='standard tab-delimited chromosome sizes file')
    parser.add_argument('output_bed', type=str, help='output TSS BED file')
    args = parser.parse_args()

    TSSCalling(**vars(args))
