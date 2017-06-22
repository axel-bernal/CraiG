#!/usr/bin/env python

import struct, string, sys, re, numpy
 
if len(sys.argv) < 3:
    sys.exit('Usage: %s coverage junctions' % sys.argv[0])
    
coverage = open(sys.argv[1], "r")
junctions = open(sys.argv[2], "r")

suffixes = ("CLmax", "CLpos", "Ymax", "Smax", "Smean", "Dscore", "HMMProb", "SigP")

if len(sys.argv) < 4 or sys.argv[3] != "add":
    for j in range(8):
        sigpout = open(prefix+"."+suffixes[j], "w")
        sigpout.close()

lengths = {}
for line in sys.stdin:
    tokens = line.split()
    lengths[tokens[1]] = int(tokens[0])

fwd_scores = {}
comp_scores = {}

for line in signalp:
    match = re.match(r"^\s*(\#|$)", line)
    if match:
        continue

    tokens = line.split()
    if len(tokens) < 15:
        print "error in format for \"", line, "\""
    
    match = re.match(r"^\s*(\S+)\.([^\.]+)\.([^\.]+)\s*$", tokens[14])

    if not match:
        print "error in format for \"", tokens[14], "\""
        break
    id = match.group(1)
    strand = match.group(2)[0:1]
    position = int(match.group(3))

    if id not in lengths:
        continue

    if id not in fwd_scores:
        fwd_scores[id] = numpy.zeros((lengths[id], 8))
        comp_scores[id] = numpy.zeros((lengths[id], 8))

    sigp = 0
    if tokens[13] == "Y":
        sigp += 1
    if tokens[20] == "Y":
        sigp += 2

    entry = (tokens[1], tokens[2], tokens[4], tokens[7], tokens[10], tokens[12], tokens[19], sigp)
    
    if strand == "0":
        fwd_scores[id][position] = entry
    else:
        comp_scores[id][lengths[id] - position + 1] = entry

for j in range(8):
    sigpout = open(prefix+"."+suffixes[j], "a")

    for id in fwd_scores:
	sigpout.write(">"+id+"\t"+str(lengths[id])+"\n")

	for i in range(lengths[id]):
	    entry = fwd_scores[id][i]
#	    print entry,"\n";
	    if entry[j] != 0:
		sigpout.write(str(i)+"\t"+str(entry[j])+"\n")

        sigpout.write("//\n")

	for i in range(lengths[id]):
	    entry = comp_scores[id][i]
#	    print entry,"\n";
	    if entry[j] != 0:
		sigpout.write(str(i)+"\t"+str(entry[j])+"\n")

    sigpout.close()
