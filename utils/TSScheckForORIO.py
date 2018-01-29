#!/usr/bin/env python

import argparse
import sys


def readRangesFromBed(args):
    def isHeader(line):
        return line.split(' ')[0] == 'track'

    ranges = []

    with open(args.input_bed) as f:
        for i, line in enumerate(f):
            if not isHeader(line):

                chromosome = line.split('\t')[0]
                start = int(line.split('\t')[1]) + 1
                end = int(line.split('\t')[2])

                if args.anchor == 'center':
                    anchor = int((start + end) / 2)
                elif args.anchor == 'start':
                    anchor = start
                elif args.anchor == 'end':
                    anchor = end

                ranges.append({
                    'bed_entry': line,
                    'chromosome': chromosome,
                    'start': anchor + args.window_start,
                    'end': anchor + args.window_start +
                    args.bin_number * args.bin_size,
                    'bed_line': i + 1,
                })

    return ranges


def readChromSizes(args):

    sizes = dict()

    with open(args.chrom_sizes) as f:
        for line in f:
            chromosome, size = line.strip().split('\t')
            sizes[chromosome] = int(size)

    return sizes


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--anchor', '-a', default='center', type=str,
                        choices=['start', 'end', 'center'], help='')
    parser.add_argument('--window_start', '-s', default=-2500, type=int,
                        help='')
    parser.add_argument('--bin_number', '-n', default=50, type=int,
                        help='')
    parser.add_argument('--bin_size', '-z', default=100, type=int,
                        help='')
    parser.add_argument('input_bed', type=str,
                        help='')
    parser.add_argument('chrom_sizes', type=str,
                        help='')
    parser.add_argument('output_bed', type=str, help='')
    args = parser.parse_args()

    input_ranges = readRangesFromBed(args)
    chrom_sizes = readChromSizes(args)

    with open(args.output_bed, 'w') as OUTPUT:
        for i, _range in enumerate(input_ranges):
            if _range['start'] > 0 and \
                    _range['end'] <= chrom_sizes[_range['chromosome']]:
                OUTPUT.write(_range['bed_entry'])
            else:
                sys.stdout.write(
                    'Entry at line ' + str(_range['bed_line']) +
                    ' runs off the chromosome. Entry removed.\n'
                )
