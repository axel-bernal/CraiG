#!/usr/bin/env python

import struct, argparse, string, sys, numpy, os, warnings, subprocess, os.path, shutil, logging, re
sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils, ioutils

parser = argparse.ArgumentParser(description='Postprocessing script for CRAIG v.1.1, a discriminative gene prediction tool for eukarya.\nWritten by Axel E. Bernal (abernal@seas.upenn.edu).')
parser.add_argument("ANNOT_FILE", type=str, help='Name of file with gene annotations', action='store')
parser.add_argument("FASTA_FILE", type=str, help='Name of file with fasta sequences or contigs', action='store')
parser.add_argument("--out-dir", type=str, help='Output directory for storing all output files.', default=".", dest='out_dir')
parser.add_argument("--use-model-scores",action='store_true', help='Uses gene model scores to break ties between overlapping gene models between samples.  By default UTR lengths are used to break ties', dest='use_gmscores');
parser.add_argument("--list-prefixes", type=str, help='Comma separated list of prefixes for the predictions.', default=".", dest='list_prefixes')
parser.add_argument("--file-prefixes", type=str, help='File with a list of prefixes for the predictions.', default=".", dest='file_prefixes')

class CraigUndefined(Exception):
    pass

try:
# checking for sanity
    if 'CRAIG_HOME' not in os.environ:
        raise CraigUndefined("CRAIG_HOME is not defined. Edit your bashrc file")

    my_env = os.environ
    my_env["PATH"] = os.environ['CRAIG_HOME']+"/perl/bin:" + \
        os.environ['CRAIG_HOME']+"/python/bin:" + \
        os.environ['CRAIG_HOME']+"/bin:" + my_env.get("PATH", '')
    
    args = parser.parse_args(sys.argv[1:])
    
    if os.path.exists(args.out_dir):
        shutil.rmtree(args.out_dir)
        
    os.makedirs(args.out_dir)

    log_file = open(args.out_dir+"/"+"postprocessing.error.log", "a")
    logging.basicConfig(filename=args.out_dir+"/"+"postprocessing.master.log", level=logging.DEBUG)
    files = []
    if args.list_prefixes != ".":
        files = args.list_prefixes.split(",") 
    elif args.file_prefixes == ".":
        raise ioutils.CraigUndefined("Must provide a valid list of prefixes") 
    else:
        files = open(args.file_prefixes,  "r").read().splitlines() 

    denovo_files = files[:];
    for i in range(len(denovo_files)):
        denovo_files[i] = denovo_files[i]+".denovo.chr"
        denovo_files[i] += ".gms" if args.use_gmscores else ""
        denovo_files[i] += ".locs"

    ioutils.clusterAnnotations(args.out_dir+"/denovo.cds", args.ANNOT_FILE, denovo_files, args.use_gmscores, my_env, log_file)

    utronly_files = files[:];
    for i in range(len(utronly_files)):
        utronly_files[i] = utronly_files[i]+".utr_only.chr"
        utronly_files[i] += ".gms" if args.use_gmscores else ""
        utronly_files[i] += ".locs"

    ioutils.clusterAnnotations(args.out_dir+"/utr_only.cds", args.out_dir+"/denovo.cds.nomatching.locs", utronly_files, args.use_gmscores, my_env, log_file)

    shell_cmd = "cat "+args.out_dir+"/denovo.cds.matching.locs "+args.out_dir+"/utr_only.cds.matching.locs "+args.out_dir+"/utr_only.cds.nomatching.locs > "+args.out_dir+"/unionized.cds.chr.locs"
    ioutils.log_and_exec(shell_cmd, my_env, log_file)

    shell_cmd = "cat "+args.out_dir+"/utr_only.cds.nomatching.locs | locs2GTF.pl "+args.FASTA_FILE+" | cat "+args.out_dir+"/denovo.cds.matching.gtf "+args.out_dir+"/utr_only.cds.matching.gtf - | gtf2gff3.pl | sed 's/5UTR/UTR/g' | sed 's/3UTR/UTR/g' | sed 's/CRAIG/CRAIG_utr_only/g' > "+args.out_dir+"/unionized.final.gff3"
    ioutils.log_and_exec(shell_cmd, my_env, log_file)

    ioutils.clusterAnnotations(args.out_dir+"/denovo", "", denovo_files, args.use_gmscores, my_env, log_file)
    shell_cmd = "cat "+args.out_dir+"/denovo.matching.gtf | gtf2gff3.pl | sed 's/5UTR/UTR/g' | sed 's/3UTR/UTR/g' | sed 's/CRAIG/CRAIG_utr_only/g' > "+args.out_dir+"/denovo.final.gff3"
    ioutils.log_and_exec(shell_cmd, my_env, log_file)


except CraigUndefined, msg:
    print msg
