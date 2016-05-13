#!/usr/bin/env python

# CREATED BY CHRISTOPHER LAVENDER
# INTEGRATIVE BIOINFORMATICS, NIEHS

import sys
import os
import argparse
from TSScall import sortList, isLessThan, isWithin


def readInBed(bed_fn, window):
    bed_list = []
    with open(bed_fn) as f:
        for line in f:
            if not ('track' in line or 'browser' in line or line == '\n'):
                chromosome, start, end, name, score, strand =\
                    line.strip().split()[:6]
                bed_list.append({
                    'chromosome': chromosome,
                    'start': int(start) - window,
                    'end': int(end) + window,
                    'name': name,
                    'score': score,
                    'strand': strand
                })
    return bed_list


def compareTo(fn_1, fn_2, list_1, list_2, missing_fn):

    if missing_fn:
        OUTPUT = open(missing_fn, 'w')
    sys.stdout.write('Bed file ' + fn_1 + ' contains ' +
                     str(len(list_1)) + ' entries.\n')

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
            if missing_fn:
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
    if missing_fn:
        OUTPUT.close()

parser = argparse.ArgumentParser()
parser.add_argument('--distance_threshold', '-d', type=int, default=0,
                    help='maximum distance threshold to be considered shared \
                    (Default: 0)')
parser.add_argument('bed_1', type=str, help='BED file 1')
parser.add_argument('bed_2', type=str, help='BED file 2')
parser.add_argument('--missing_entries_file', type=str, default=None,
                    help='output file with missing entries')
args = parser.parse_args()

assert os.path.exists(args.bed_1)
assert os.path.exists(args.bed_2)

bed_list_1 = readInBed(args.bed_1, 0)
bed_list_2 = readInBed(args.bed_2, args.distance_threshold)

compareTo(args.bed_1, args.bed_2, bed_list_1, bed_list_2,
          args.missing_entries_file)
