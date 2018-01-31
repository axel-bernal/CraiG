#!/usr/bin/env python
import sys, argparse, re, warnings, numpy
import tempfile, subprocess, os
import traceback
import csv

parser = argparse.ArgumentParser(description='Convert coverage into score filter. Coverage files follow UCSC format for coordinates.')
parser.add_argument('INPUT_FILE', type=str, help='input file', action='store')
parser.add_argument('--len-file', type=argparse.FileType('r'),
                    help='file with contig lengths', dest='len_file')
parser.add_argument('--format', type=str, help='format for input file',
                    choices=['rum', 'bam'], default='bam')
parser.add_argument('-o', '--orientation', type=str, help='RNA-Seq reads are stranded and thus have orientation. Possible values are N,F,R,RF,FR, FF or NN. A value of N indicates no orientation and thus no strand can be extracted. Only FR (Illumina) supported for now', dest='orientation', default='FR')

line = ""
if len(sys.argv) < 2:
    sys.exit('Usage: %s --format FMT --len-file LEN_FILE < coverage' % sys.argv[0])

lens = {}
coverages = {}
min_coverage = 1.0
record = numpy.dtype([('cov', numpy.float64), ('end', numpy.int32)])

def process_coverage(fh, strand):
    global coverages, lens, min_coverage
    last_start = last_contig = last_coverage = None
    blksize = 1
    for tokens in csv.reader(fh, delimiter="\t"):
        if args.format == 'rum':
            coverage = float(tokens[3])
        else:
            coverage = float(tokens[2])

        if coverage == 0:
            continue
        if coverage < min_coverage:
            min_coverage = coverage

        contig = tokens[0]
        if contig not in lens:
            raise Exception("length for contig id {} NOT FOUND".format(contig))
        if contig not in coverages:
            coverages[contig] = (numpy.zeros(lens[contig] + 1, dtype=record),
                                 numpy.zeros(lens[contig] + 1, dtype=record))

        if args.format == 'rum':
            start = int(tokens[1]) + 1;
        else:
            start = int(tokens[1])
            if not last_start is None and last_start == start - blksize and \
               last_contig == contig and last_coverage == coverage:
                start = last_start
                blksize += 1
            else:
                blksize = 1
                if last_contig is not None and contig != last_contig:
                    last_start = None
                else:
                    last_start = start
                last_contig = contig
                last_coverage = coverage

        if args.format == 'rum':
            end = int(tokens[2])
        else:
            end = int(tokens[1])

        coverages[contig][strand][start]['cov'] = coverage
        coverages[contig][strand][start]['end'] = end


try:
    args = parser.parse_args(sys.argv[1:])
    inpfiles = [args.INPUT_FILE]
    non_stranded = args.orientation == "N" or args.orientation == 'NN'
    script_dir = os.path.abspath(os.path.dirname(__file__))

    for line in args.len_file:
        match = re.match(r"^(\d+)\s+(\S+)$", line)
        if match:
            lens[match.group(2)] = int(match.group(1))

    if args.format == 'bam':
        if args.INPUT_FILE == '-':
            raise Exception("Cannot stream BAMs from stdin. Must provide name")

        inpfiles[0] = (tempfile.NamedTemporaryFile(delete=False)).name
        if args.orientation == 'FR':
            inpfiles.append((tempfile.NamedTemporaryFile(delete=False)).name)
            f2reads = 'samtools view -u -f 128 -F 16 {} | samtools depth /dev/stdin | grep -v "\\t0$"'.\
                      format(args.INPUT_FILE)
            r1reads = 'samtools view -u -f 80 {} | samtools depth /dev/stdin | grep -v "\\t0$"'.\
                      format(args.INPUT_FILE)
            merge_script = "%s/merge_samtools_depths.py --len-file %s" % \
                           (script_dir, args.len_file.name)
            shell_cmd = '{ %s ; %s ; } | %s > %s' % \
                        (f2reads, r1reads, merge_script, inpfiles[0])
            sys.stderr.write("Running " + shell_cmd + "\n")
            subprocess.check_call(shell_cmd, shell=True)
            r2reads = 'samtools view -u -f 144 {} | samtools depth /dev/stdin | grep -v "\\t0$"'.\
                      format(args.INPUT_FILE)
            f1reads = 'samtools view -u -f 64 -F 16 {} | samtools depth /dev/stdin | grep -v "\\t0$"'.\
                      format(args.INPUT_FILE)
            shell_cmd = '{ %s ; %s ; } | %s > %s' % (r2reads, f1reads, merge_script, inpfiles[1])
            sys.stderr.write("Running " + shell_cmd + "\n")
            #subprocess.check_call(shell_cmd, shell=True)
        elif non_stranded:
            shell_cmd = "samtools depth {} > {}".format(args.INPUT_FILE, inpfiles[0])
            subprocess.check_call(shell_cmd, shell=True)
        else:
            raise Exception("Orientation {} unavailable".format(args.orientation))
    elif not non_stranded:
        raise Exception("Cannot process stranded RUM alignments")

    for strand in range(len(inpfiles)):
        with open(inpfiles[strand], "r") if inpfiles[strand] != '-' else sys.stdin as fh:
            process_coverage(fh, strand)

    for contig in coverages:
        sys.stdout.write(">"+contig+"\t"+str(lens[contig])+"\n")
        for strand in range(0,2):
            for i in range(1,lens[contig] + 1):
                coverage = coverages[contig][strand][i]['cov']
                if coverage == 0:
                    continue

                coverage = int(coverage/min_coverage)
                end = coverages[contig][strand][i]['end']
                pend = end
                pstart = i
                if args.format == 'bam' and strand:
                    pstart = lens[contig] - end + 1
                    pend = lens[contig] - i + 1

                buffer = str(pstart)
                if i < end:
                    buffer += ".." + str(pend)
                buffer += "\t" + str(coverage) + "\n"
                sys.stdout.write(buffer)

            if strand == 0 and not non_stranded:
                sys.stdout.write("//\n")
            else:
                break

except Exception, e:
    print "Exception in user code:"
    print '-'*60
    traceback.print_exc(file=sys.stderr)
    print '-'*60
