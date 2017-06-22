#!/usr/bin/env python

import struct, argparse, string, sys, numpy, os
sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils
 
if len(sys.argv) < 2:
    sys.exit('Usage: %s FILE [-by <i|c|l>] [-sep RECORD_SEP]' % sys.argv[0])
    
line = ""
parser = argparse.ArgumentParser(description='sorts list of records using id, contig or location as key.')
parser.add_argument("FILE", type=str, help='First input filename', action='store')
parser.add_argument("-by", type=str, help='<i|c|l>', default='i')
parser.add_argument("-sep", type=str, help='record separator', default="\\n>")

records = {}
idType = 'i'

try:
    args = parser.parse_args(sys.argv[1:])
    idType = args.by

    fh = open(args.FILE, "r") if args.FILE != '-' else sys.stdin

    for line in bioutils.readRecords(fh, args.sep.decode("string_escape")):
        (key, value) = bioutils.processRecord(line, idType)
        if key != None:
            if key in records:
                records[key].append(value);
            else:
                records[key] = [value]
        
    if args.FILE != '-': fh.close()

except AttributeError: sys.stdout.write("Error reading \""+line+"\"")

try:

    for key in sorted(records.iterkeys()):
        for value in records[key]:
            bioutils.printRecord(key, value, args.sep, idType)
except AttributeError: sys.stdout.write("Error writing \""+key+"\"")

