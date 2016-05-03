#!/usr/bin/env python

# CREATED BY CHRISTOPHER LAVENDER
# INTEGRATIVE BIOINFORMATICS, NIEHS

import sys
import os
import argparse
from TSScall import sortList, isLessThan, isWithin


def readInBed(bed_fn):
    bed_list = []
    with open(bed_fn) as f:
        for line in f:
            if not ('track' in line or 'browser' in line or line == '\n'):
                chromosome, start, end, name, score, strand =\
                    line.strip().split()[:6]
                bed_list.append({
                    'chromosome': chromosome,
                    'start': int(start),
                    'end': int(end),
                    'name': name,
                    'score': score,
                    'strand': strand
                })
    return bed_list


def compareTo(fn_1, fn_2, list_1, list_2, missing_fn):
    with open(missing_fn, 'w') as OUTPUT:
        sys.stdout.write('Bed file ' + fn_1 + ' contains ' +
                         str(len(list_1)) + ' entries.\n')

        # position_list_1 = []
        # list_2_positions = []
        # shared = 0
        #
        # names_and_positions_list_1 = []
        # for entry in list_1:
        #     names_and_positions_list_1.append({
        #         'name': entry['name'],
        #         'position':entry['strand']+':'+entry['chromosome']+':'+str(entry['start'])+':'+str(entry['end'])
        #         })
        #
        # positions_list_2 = []
        # for entry in list_2:
        #     positions_list_2.append(entry['strand']+':'+entry['chromosome']+':'+str(entry['start'])+':'+str(entry['end']))
        #
        # for entry in names_and_positions_list_1:
        #     if entry['position'] in positions_list_2:
        #         shared += 1
        #     else:
        #         OUTPUT.write(entry['name'] + '\t' + entry['position'] + '\n')

        list_1 = sortList(list_1, 'sort_by_strand')
        list_2 = sortList(list_2, 'sort_by_strand')

        shared = 0

        list_1_index = 0
        list_2_index = 0

        while list_1_index < len(list_1) and list_2_index < len(list_2):
            if isWithin(list_1[list_1_index], list_2[list_2_index]):
                shared += 1
                list_1_index += 1
            elif isLessThan(list_1[list_1_index], list_2[list_2_index]):
                OUTPUT.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    list_1[list_1_index]['chromosome'],
                    list_1[list_1_index]['start'],
                    list_1[list_1_index]['end'],
                    list_1[list_1_index]['name'],
                    list_1[list_1_index]['score'],
                    list_1[list_1_index]['strand']
                    ))
                list_1_index += 1
            else:
                list_2_index += 1

    sys.stdout.write(str(shared) + ' are shared with ' + fn_2 + '.\n')

parser = argparse.ArgumentParser()
parser.add_argument('bed_1', type=str, help='BED file 1')
parser.add_argument('bed_2', type=str, help='BED file 2')
parser.add_argument('missing_entries_file', type=str,
                    help='output file with missing entries')
args = parser.parse_args()

assert os.path.exists(args.bed_1)
assert os.path.exists(args.bed_2)

bed_list_1 = readInBed(args.bed_1)
bed_list_2 = readInBed(args.bed_2)

compareTo(args.bed_1, args.bed_2, bed_list_1, bed_list_2,
          args.missing_entries_file)
