#!/usr/bin/env python
import sys
import argparse
import csv
import re
import numpy

parser = argparse.ArgumentParser(description='Merges two samtools depth outputs to stdout')
#parser.add_argument('DEPTH_FILE1', type=str, help='First depth file', action='store')
#parser.add_argument('DEPTH_FILE2', type=str, help='Second depth file', action='store')
parser.add_argument('--len-file', type=argparse.FileType('r'),
                    help='file with contig lengths', dest='len_file')

depths = {}
lens = {}

# Aggregate sales by contig, position
def aggregate_depth(file):
    global depths, lens
    for row in csv.reader(file, delimiter="\t"):
        contig = row[0]
        if contig not in lens:
            raise Exception("length for contig id {} NOT FOUND".format(contig))
        if contig not in depths:
            depths[contig] = (numpy.zeros(lens[contig] + 1, dtype=int))

        depths[contig][int(row[1])] += int(row[2])


args = parser.parse_args(sys.argv[1:])

for line in args.len_file:
    match = re.match(r"^(\d+)\s+(\S+)$", line)
    if match:
        lens[match.group(2)] = int(match.group(1))

#read_depth(open(args.DEPTH_FILE1)
#read_depth(open(args.DEPTH_FILE2))
aggregate_depth(sys.stdin)

# Write aggregated data to file.

for contig in depths:
    for pos in xrange(lens[contig] + 1):
        depth = depths[contig][pos]
        if not depth:
            continue
        print "{}\t{}\t{}".format(contig, pos, depth)
