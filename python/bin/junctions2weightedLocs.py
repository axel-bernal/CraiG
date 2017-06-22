#!/usr/bin/env python
import sys
import argparse
import re
import warnings
import tempfile
import subprocess

parser = argparse.ArgumentParser(description='Convert junctions into weighted locations.')
parser.add_argument('INPUT_FILE', type=str, help='input file', action='store')
parser.add_argument('OUTPUT_FILE', type=argparse.FileType('w'), help='output file')
parser.add_argument('--dataset-id', type=str,
                    help='Dataset description. Make sre it is spelled correctly as it will be used for matching against the corresponding fields in case of epdb* format', dest='dataset_id', default="dataset1")
parser.add_argument('--sample-id', type=str,
                    help='Sample description. See --dataset-id option for details',
                    dest='sample_id', default="")
parser.add_argument('--with-bed', type=str,
                    help='Generates bed file with junctions', default="")

parser.add_argument('--cut-off', type=int,
                    help='RNA junction coverage cut-off', dest='score_cutoff', default=2)
parser.add_argument('--format', type=str, choices=['bam', 'rum', 'gsnap', 'epdb1', 'epdb2'],
                    help='one of epdb1(plasmo), epdb2(toxo), rum, bam or gsnap',
                    default='bam')
parser.add_argument('--min-intron-len', type=int,
                    help="Minimum intron length", default=20)
parser.add_argument('-o', '--orientation', type=str,
                    help='RNA-Seq reads are stranded and thus have orientation. Possible values are N,F,R,RF,FR, FF or NN. A value of N indicates no orientation and thus no strand can be extracted', dest='orientation', default='N')


def mapping_coords(coords):
    pos = coords.find(":")
    mapping = coords[pos + 1:].split("..")
    mapping.insert(0, coords[:pos])
    return mapping


try:
    args = parser.parse_args(sys.argv[1:])
    inpfile = args.INPUT_FILE
    if args.format == 'bam':
        if inpfile == '-':
            raise Exception("Cannot stream BAMs from stdin. Must provide name")
        if args.with_bed:
            bedfile = args.with_bed
        else:
            bedfile = (tempfile.NamedTemporaryFile(delete=False)).name
        shell_cmd = "regtools junctions extract -i {} {} > {}" \
                    .format(args.min_intron_len, inpfile, bedfile)

        subprocess.check_call(shell_cmd, shell=True)
        inpfile = bedfile
    elif args.with_bed:
        raise Exception("Bed file output not implemented for {} format".format(args.format))

    counter = 1
    state = 0
    read_no = 0
    nsegs = 0

    with open(inpfile) as f_:
        for line in f_:
            if args.format != 'gsnap':
                match = re.match(r"^\s*\#", line)
                if match:
                    continue

                tokens = line.split()
                if args.format == 'rum':
                    score = int(tokens[1])
                elif args.format == 'bam':
                    score = int(tokens[4])
                else:
                    score = int(tokens[6])

                if score < args.score_cutoff:
                    continue

                contig = tokens[0] if args.format in ['epdb2', 'bam'] else tokens[3]
                if args.format == 'bam':
                    blksizes = tokens[10].split(",")
                    strand = tokens[5]
                    if strand == '+':
                        start = str(int(tokens[1]) + int(blksizes[0]))
                        end = str(int(tokens[2]) - int(blksizes[1]) + 1)
                    else:
                        start = str(int(tokens[2]) - int(blksizes[1]) + 1)
                        end = str(int(tokens[1]) + int(blksizes[0]))
                else:
                    start = str(int(tokens[4]) - 1)
                    end = str(int(tokens[5]) + 1)

                if args.format == 'bam':
                    ident = tokens[3]
                else:
                    if args.format == 'epdb2':
                        ident = tokens[3]+"."+str(counter)
                    else:
                        ident = tokens[2]+"."+tokens[0]
                    ident += "."+tokens[7]

                if args.format not in ['epdb1', 'epdb2', 'bam']:
                    match = re.match(r"^([^\:]+)\:(\d+)\-(\d+)", tokens[0])
                    contig = match.group(1)
                    start = str(int(match.group(2)) - 1)
                    end = str(int(match.group(3)) + 1)
                    ident = args.dataset_id+args.sample_id+"."+str(counter)+"."+tokens[2]

                args.OUTPUT_FILE.write(">"+ident+"\t<"+contig+"_"+start+"_"+start+";")
                args.OUTPUT_FILE.write(contig+"_"+end+"_"+end+">\t"+str(score)+";"+str(score)+"\n")
            else:
    #            print state, len(line), line
                if not line.strip():
                    continue

                tokens = line.split()
    #            print >> sys.stderr, "READ ", read_no, state, nsegs, line
    #            print nsegs, state, line
                if tokens[0][0] == ">" or tokens[0][0] == "<":
                    state = 1
                    read_no = (0 if tokens[0][0] == ">" else 1)

                    if nsegs != 0:
                        sys.exit("number of segments should be zero but it's "+str(nsegs)+" instead")
                    if read_no and len(args.orientation) < 2:
                        sys.exit("orientation should be given for each paired read "+args.orientation)
                    continue

                if state == 1 and len(tokens) >= 5:
                    details = tokens[3].split(",")
                    if len(details) < 4:
                        continue

                    direction = details[3]
                    match = (direction=="dir:antisense" or direction=="dir:sense")
    #                match = re.match(r"^ \S+\s+\S+\s+(\-|\+)?([^\:]+)\:(\d+)\.\.(\d+).+(acceptor|donor).+dir\:(sense|antisense).+\s+segs\:(\d+)", line)

                    if match:
                        cmap = mapping_coords(tokens[2])
                        cstrand = cmap[0][0]
                        contig = cmap[0][1:]
                        beg = int(cmap[1])
                        end = int(cmap[2])
                        state = 2
                        nsegs = int(tokens[4][5]) - 1
                    continue

                if state == 2 and nsegs >= 1:
                    prevcontig = contig
                    prevcstrand = cstrand
                    prevbeg = beg
                    prevend = end
                    prevdirection = direction
    #                print >> sys.stderr, "decreasing nsegs ", nsegs
                    nsegs -= 1
                    if nsegs == 0:
                        state = 1

                    if len(tokens) < 4:
                        continue

                    details = tokens[3].split(",")
                    if len(details) < 4:
                        continue

                    direction = details[3]
                    match = (direction=="dir:antisense" or direction=="dir:sense")
    #                match = re.match(r"^,\S+\s+\S+\s+(\-|\+)?([^\:]+)\:(\d+)\.\.(\d+).+(acceptor|donor).+dir\:(sense|antisense)", line)
                    if match:
                        cmap = mapping_coords(tokens[2])
                        cstrand = cmap[0][0]
                        contig = cmap[0][1:]
                        beg = int(cmap[1])
                        end = int(cmap[2])

                        if prevdirection != direction:
                            raise Exception, "direction of aligned segments for contig "+contig+" are different "+prevdirection+" and "+direction
                        if prevcstrand != cstrand or prevcontig != contig :
                            warnings.warn("strands or contigs of aligned segments "+prevcontig+" and "+contig+" should be the same but are different. At line "+line)
                            warnings.warn("A potential translocation splicing event")
                            prevcontig = contig
                            prevcstrand = cstrand
                            prevbeg = beg
                            prevend = end
                            prevdirection = direction
                            continue

                        ident = args.dataset_id+args.sample_id+"."+str(counter)

                        ibeg = prevend
                        iend = beg
                        offset = 1
                        if direction == "dir:antisense":
                            iend = prevend
                            ibeg = beg
                            if args.orientation[read_no] == "F":
                                warnings.warn("splice orientation for "+contig +"_"+str(ibeg)+"_"+str(iend)+" is different to what the given orientation suggest "+args.orientation[read_no]+" "+direction)
                                continue
                        else:
                            if args.orientation[read_no] == "R":
                                warnings.warn("splice orientation for "+contig+"_"+str(ibeg)+"_"+str(iend)+" is different to what the given orientation suggest "+args.orientation[read_no]+" "+direction)
                                continue
                        if ibeg > iend:
                            offset = -1

                        args.OUTPUT_FILE.write(">"+ident+"\t<"+contig + "_" + str(ibeg - offset) + "_" + str(ibeg) + ";")
                        args.OUTPUT_FILE.write(contig + "_" + str(iend) + "_" +
                                         str(iend + offset) + ">\t" + "1;1\n")


            counter += 1

        if args.OUTPUT_FILE.name != "-":
            args.OUTPUT_FILE.close()

except Exception, e:
    sys.stderr.write(str(e))
