#!/usr/bin/env python

import argparse
import sys
import os

sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils

parser = argparse.ArgumentParser(description='intersects 2 lists of records using id, contig or location as key.')
parser.add_argument("FILE1", type=str, help='First input filename',
                    action='store')
parser.add_argument("FILE2", type=str, help='Second input filename',
                    action='store')
parser.add_argument("-by", type=str, help='A combination of <i|c|l><i|c|l>',
                    default='ii')
parser.add_argument("-sep", type=str, help='record separator', default="\\n>")
parser.add_argument("--num-fields1", type=int,
                    help='number of fields(columns) for records in FILE1',
                    default=1)
parser.add_argument("--num-fields2", type=int,
                    help='number of fields(columns) for records in FILE2',
                    default=1)
parser.add_argument("-v", type=int, help='suppress output columns from FILEV',
                    default=2)
parser.add_argument("-u", "--print-unpaired", action='store_true',
                    help='print unpaired lines with blank field',
                    dest='print_unpaired')

line = ""
records = {}
idType = ['i', 'i']
empty_recval1 = ""
empty_recval2 = ""

try:
    args = parser.parse_args(sys.argv[1:])
    for i in range(2): idType[i] = args.by[i:i+1]

    for i in range(args.num_fields1):
        empty_recval1 += "\t "
    for i in range(args.num_fields2):
        empty_recval2 += "\t "

    fh = open(args.FILE1, "r") if args.FILE1 != '-' else sys.stdin

    for line in bioutils.readRecords(fh, args.sep.decode("string_escape")):
        (key, value) = bioutils.processRecord(line, idType[0])
        if key == None:
            continue

        if key in records: records[key][0].append(value)
        else:  records[key] = [[value], [empty_recval1]]

    if args.FILE1 != '-': fh.close()

    fh = open(args.FILE2, "r") if args.FILE2 != '-' else sys.stdin

    for line in bioutils.readRecords(fh, args.sep.decode("string_escape")):
        (key, value) = bioutils.processRecord(line, idType[1])
        if key == None:
            continue

        if key in records:
            if records[key][1][0] == empty_recval1:
                records[key][1][0] = value
            else:
                records[key][1].append(value)
        else: records[key] = [[empty_recval2], [value]]

    if args.FILE2 != '-': fh.close()

    for key in records:
        unpaired = (records[key][0][0] == empty_recval1 or records[key][1][0] == empty_recval2)
        if not args.print_unpaired and unpaired:
            continue

        for i in range(len(records[key][0])):
            final_value = records[key][0][i] if args.v != 1 else ""
#            print "\"AR1=", final_value, "\""
            for j in range(len(records[key][1])):
#                print "Record \"", records[key][0][i], records[key][1][j], "\""
                final_value += records[key][1][j] if args.v != 2 else ""
#                print "\"FV=", final_value, "\"","\n"
                bioutils.printRecord(key, final_value, args.sep, idType[0])

except (AttributeError, TypeError, RuntimeError, NameError): print line, " did not match record\n"
