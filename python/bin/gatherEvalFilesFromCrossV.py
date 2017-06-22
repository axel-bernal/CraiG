#!/usr/bin/env python
# input to the program are head -15 of eval files. Something like this
# ls cross1.?/*.eval | grep -v "SP"| xargs head -15 $1 | gatherEvalFilesFromCrossV.py cross/test.locs
#
import struct, string, sys, re

if len(sys.argv) < 2:
    sys.exit('Usage: %s output-file' % sys.argv[0])

outfile = open(sys.argv[1], "w");
entries = {}
exp = -1;
sn = 0;
locFile = "";
score = 0;

for line in sys.stdin:

    match = re.match(r"^==>\s(\S+)\.(\d)\/([\w,\.]+)\.eval", line)
    if match:
        exp = match.group(2)
        locFile = match.group(1) + "." + exp + "/" + match.group(3) + ".locs"
        continue

    match = re.match(r"^Gene\sSensitivity\s+([^\%]+)\%", line)
    if match:
        sn = float(match.group(1))
        continue

    match = re.match(r"^Gene\sSpecificity\s+([^\%]+)\%", line)
    if match:
        sp = float(match.group(1))
        if sn+sp == 0:
            score = 0
        else:
            score = 2*sn*sp/(sn+sp)
        
#        print locFile, ' ', score, ' ', sn, ' ', sp, '\n'
        if exp in entries:
            if entries[exp][1] < score:
                entries[exp] = (locFile, score);
#                print "pick ", locFile, ' ', score, ' ', sn, ' ', sp, '\n'
        else:
            entries[exp] = (locFile, score)
#            print "pick ", locFile, ' ', score, ' ', sn, ' ', sp, '\n'
        continue

for exp in entries:
    infile = open(entries[exp][0],"r")
    outfile.writelines(infile.readlines())
    infile.close()

outfile.close()

    
