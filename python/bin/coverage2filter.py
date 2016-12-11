#!/usr/bin/env python                                                                
import sys, argparse, re, warnings, numpy

parser = argparse.ArgumentParser(description='Convert coverage into score filter. Coverage files follow UCSC format for coordinates.')
parser.add_argument('--len-file', type=argparse.FileType('r'), help='file with contig lengths', dest='len_file')
parser.add_argument('--format', type=str, help='one of epdb or rum', default='epdb')
parser.add_argument('-s', '--stranded', action='store_true', help='RNA-Seq reads provide strand information', dest='stranded')
line = ""

if len(sys.argv) < 2:
    sys.exit('Usage: %s --format FMT --len-file LEN_FILE < coverage' % sys.argv[0])

try:
    lens = {}
    coverages = {}
    args = parser.parse_args(sys.argv[1:])
    record = numpy.dtype([('cov', numpy.float64), ('end', numpy.int32)])

    for line in args.len_file:
        match = re.match(r"^(\d+)\s+(\S+)$", line)
        lens[match.group(2)] = int(match.group(1))

    min_coverage = 1.0

    for line in sys.stdin:
        match = re.match(r"^\s*\#", line)
        if match:  continue
        
        tokens = line.split('\t')

        multiple = int(tokens[8]) if args.format == 'epdb' else 0;
        if multiple == 1:
            continue

        coverage = float(tokens[7]) if args.format == 'epdb' else float(tokens[3]);
        if coverage < min_coverage:
            min_coverage = coverage

        contig = tokens[0]

        if contig not in coverages:
            coverages[contig] = (numpy.zeros(lens[contig] + 1, dtype=record),
                                 numpy.zeros(lens[contig] + 1, dtype=record))
                                 
        strand = 0
        if args.format == 'epdb':
            if tokens[6] != '':
                strand = int(tokens[6])
                args.stranded = strand

        if not args.stranded and strand != 0:
            raise Exception, "non-stranded RNASeq with complementary strand elements"

        start = int(tokens[4]) + 1 if args.format == 'epdb' else int(tokens[1]) + 1;
        end = int(tokens[5]) if args.format == 'epdb' else int(tokens[2]);

        coverages[contig][strand][start]['cov'] = coverage
        coverages[contig][strand][start]['end'] = end
        
    for contig in coverages:
        sys.stdout.write(">"+contig+"\t"+str(lens[contig])+"\n")
        for strand in range(0,2):
            for i in range(1,lens[contig] + 1):
                coverage = coverages[contig][strand][i]['cov']
                if coverage == 0:
                    continue

                coverage = int(coverage/min_coverage)
                sys.stdout.write(str(i))
                end = coverages[contig][strand][i]['end']
                if i < end:
                    sys.stdout.write(".."+str(end))
                sys.stdout.write("\t"+str(coverage)+"\n")

            if strand == 0 and args.stranded: sys.stdout.write("//\n")
            else : break

except (IndexError,AttributeError): sys.stderr.write(line+" did not match record\n") 
