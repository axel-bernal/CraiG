#!/usr/bin/env python
import sys, argparse, re, warnings, numpy

parser = argparse.ArgumentParser(description='Convert gsnap reads to coverage filter.')
parser.add_argument('--len-file', type=argparse.FileType('r'), help='file with contig lengths', dest='len_file')
parser.add_argument('-nu','--non-unique', action='store_true', help='Report non-unique mappings. By default only unique mappings are reported', dest='nunique')
parser.add_argument('-o', '--orientation', type=str, help='RNA-Seq reads are stranded and thus have orientation. Possible values are N,F,R,RF,FR, FF or NN. A value of N indicates no orientation and thus no strand can be extracted', dest='orientation', default='NN')

if len(sys.argv) < 2:
    sys.exit('Bad Usage: use '+sys.argv[0] + ' -h for options')

def mapping_coords(coords):
    pos = coords.find(":")
    mapping = coords[pos + 1:].split("..")
    mapping.insert(0, coords[:pos])
    return mapping

counter = 1
read_no = 0
try:
    lens = {}
    coverages = {}
    args = parser.parse_args(sys.argv[1:])

    for line in args.len_file:
        match = re.match(r"^(\d+)\s+(\S+)$", line)
        lens[match.group(2)] = int(match.group(1))

    min_coverage = 1.0

    for line in sys.stdin:
        if not line.strip():
            continue

        tokens = line.split()

        if tokens[0][0]==">" or tokens[0][0]=="<":
            read_no = (0 if tokens[0][0] == ">" else 1)
            if read_no and len(args.orientation) < 2:
                sys.exit("orientation should be given for each paired read")

            unique = True if tokens[1] == "1" else False
            continue

        if unique == args.nunique:
            continue

        cmap = mapping_coords(tokens[2])
        cstrand = cmap[0][0]
        contig = cmap[0][1:]
        beg = int(cmap[1])
        end = int(cmap[2])
#        print contig, beg, end
#        match = re.match(r"^(\s|,)\S+\s+\S+\s+(\-|\+)?([^\:]+)\:(\d+)\.\.(\d+)\s+\S.+", line)
#        if not match:
#            continue

        if contig not in coverages:
            coverages[contig] = (numpy.zeros(lens[contig]+1, dtype=numpy.int),
                                 numpy.zeros(lens[contig]+1, dtype=numpy.int))

        strand = 0
        if args.orientation[read_no] == "R":
            beg = end
            end = int(cmap[1])

        if beg > end:
            if args.orientation[read_no] != "N":
                strand = 1
                beg = lens[contig] - beg + 1
                end = lens[contig] - end + 1
            else:
                tmp = beg
                beg = end
                end = tmp

        for i in range(beg, end + 1):
            coverages[contig][strand][i] += 1

    for contig in coverages:
        sys.stdout.write(">"+contig+"\t"+str(lens[contig])+"\n")
        for strand in range(0,2):
            last_pos = 1
            last_cov = coverages[contig][strand][last_pos]
            for i in range(2,lens[contig] + 1):
                cov = coverages[contig][strand][i]
                if cov != last_cov:
                    if last_cov != 0:
                        sys.stdout.write(str(last_pos))
                        end = i - 1
                        if last_pos < end:
                            sys.stdout.write(".." + str(end))
                        sys.stdout.write("\t" + str(last_cov) + "\n")
                    last_cov = cov
                    last_pos = i

            if last_cov != 0:
                sys.stdout.write(str(last_pos))
                if last_pos < lens[contig]:
                    sys.stdout.write(".." + str(lens[contig]))
                sys.stdout.write("\t" + str(last_cov) + "\n")

            if strand==0 and args.orientation!='N' and args.orientation!='NN':
                sys.stdout.write("//\n")
            else : break
except AttributeError: sys.stderr.write(line+" did not match record\n")
