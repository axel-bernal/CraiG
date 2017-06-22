#!/usr/bin/env python

import struct, argparse, string, sys, numpy, os, warnings
sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils
  
parser = argparse.ArgumentParser(description='computes the set difference between 2 lists of records using id, contig or location as key.')
parser.add_argument("FILE1", type=str, help='First input filename', action='store')
parser.add_argument("FILE2", type=str, help='Second input filename', action='store')
parser.add_argument("-by", type=str, help='A combination of <i|c|l><i|c|l>', default='ii')
parser.add_argument("-sep", type=str, help='record separator', default="\\n>")

line = ""
records = {}
idType = ['i', 'i']

try:   
    args = parser.parse_args(sys.argv[1:])
    for i in range(2): idType[i] = args.by[i:i+1]

    fh = open(args.FILE2, "r") if args.FILE2 != '-' else sys.stdin

    for line in bioutils.readRecords(fh, args.sep.decode("string_escape")):
        (key, value) = bioutils.processRecord(line, idType[1])
        if key != None: records[key] = 1;

    if args.FILE2 != '-': fh.close()
    
    fh = open(args.FILE1, "r") if args.FILE1 != '-' else sys.stdin

    for line in bioutils.readRecords(fh, args.sep.decode("string_escape")):
        (key, value) = bioutils.processRecord(line, idType[0])
        if key != None and key not in records:
            bioutils.printRecord(key, value, args.sep, idType[0])
    if args.FILE1 != '-': fh.close()
except AttributeError: print line, " did not match record\n"
