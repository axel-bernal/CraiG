#!/usr/bin/env python

import struct, string, sys, numpy, os, warnings, re
sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils
line = ""
clusters = {}
sep = "\\n>"
exonParse1 = re.compile(r"\[")
exonParse2 = re.compile(r"\]")

class CraigUndefined(Exception):
    pass

try:
    for i in range(1, len(sys.argv)):
        filen = sys.argv[i]
        fh = open(filen, "r") if filen != '-' else sys.stdin

        for line in bioutils.readRecords(fh, sep.decode("string_escape")):
            value = line.rstrip('\n')
            tokens = value.split()
            fexons = exons = lexons = []
            bioutils.readLocation(tokens[1], fexons, exons, lexons)
            print fexons, exons, lexons
            if len(exons) == 0:
                continue

            if len(fexons) > 1 or len(exons) > 1 or len(lexons) > 1:
                if len(lexons) > 0:
                    exons[-1][2] = lexons[0][2]
                    for i in range(1, len(lexons)):
                        exons.append(lexons[i])
                if len(fexons) > 0:
                    exons[0][1] = fexons[-1][1]
                    for i in range(len(fexons) - 2, 0, -1):
                        exons.insert(0, fexons[i])

                exons[0][1] = exons[0][2]
                exons[-1][2] = exons[0][1]
                key = "_".join(exons[0])
                for i in range(1, len(exons)):
                    key = key + ";" + "_".join(exons[i])
            else:
                key = exons[0][0]+"_"+exons[0][2]

            print key, "->", value, "\n";

            if key != None:
                base_filen = os.path.basename(filen)
                if key in clusters:
                    clusters[key].extend([value, filen]);
                else:
                    clusters[key] = [value, filen]

        if filen != '-': fh.close()
    
    for key, locs in clusters:
        for i in range(0, len(locs), 2):
            print ">" + locs[i] + " " + locs[i+1] + "\n"
        print "//\n"

except CraigUndefined, msg:
    print msg
