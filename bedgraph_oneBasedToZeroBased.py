#!/usr/bin/env python

import sys

def printUsage():
    sys.stdout.write("Usage: python bedgraph_zeroBased2oneBased.py <input_file> <output_file>\n")
    sys.exit(2)

if len(sys.argv) != 3:
    printUsage()
    
OUTPUT = open(sys.argv[2], "w")

with open(sys.argv[1]) as f:
    for line in f:
        #if line != "":
            if line[0] != "#":
                line_split = line.strip().split()
                
                if len(line_split) != 0:
                    if line_split[0] == "track":
                        OUTPUT.write(line + "\n")
                    elif line_split[0] != "track":
                        OUTPUT.write(line_split[0] + "\t" + str(int(line_split[1])-1) + "\t" + line_split[2] + "\t" + line_split[3] + "\n")