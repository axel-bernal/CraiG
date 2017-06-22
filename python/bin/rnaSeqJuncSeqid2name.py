#!/usr/bin/env python                                                                
import struct, string, sys, re, numpy

idmapfile = open(sys.argv[1], "r")

seqids = {}
for line in idmapfile:
    tokens = line.split()
    seqids[tokens[0]] = tokens[1]

for line in sys.stdin:
    tokens = line.split()
    tokens[3] = seqids[tokens[3]]
    print "\t".join(tokens)


