#!/usr/bin/env python

import struct, argparse, string, sys, numpy, os, warnings, subprocess, os.path, shutil, logging, re
sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils, ioutils
parser = argparse.ArgumentParser(description='Script for building task.prop, conttroller.prop and all the associated directories for executing the CRAIG v.1.1 pipeline, a discriminative gene prediction tool for eukarya.\nWritten by Axel E. Bernal (abernal@seas.upenn.edu).')

parser.add_argument("SPECIES", type=str, help='Species name, given by the common suffix of files present in the $CRAIG_HOME/models directory. These files specify the resources, filters and features used for constructing gene models and will be generated. Possible values right now are "tgondii" "cvelia" and "vitrella"', action='store')
parser.add_argument("ANNOT_FILE", type=str, help='Name of file with gene annotations', action='store')
parser.add_argument("FASTA_FILE", type=str, help='Name of file with fasta sequences or contigs', action='store')
parser.add_argument("RSQ_INPUTDIR", type=str, help='Mapped RNA-Seq input directory, i.e. where the RNA-Seq mapped data is located"', action='store')
parser.add_argument("PREFIX_CRAIG_OUTPUT", type=str, help='Prefix for CRAIG output files. Output will be found in PREFIX_CRAIG_OUTPUT.preproc and PREFIX_CRAIG_OUTPUT.model. Temporary files for DJob execution will be found in PREFIX_CRAIG_OUTPUT.test. The latter should be erased manually after results are obtained', action='store')
parser.add_argument("PREFIX_PROPDIR", type=str, help='Distribjob input directory, i.e. where task.prop and controller.prop are going to be located', action='store')
parser.add_argument("--transcript-tag", type=str, help='The description tag as it appears in ANNOT_FILE for transcript regions(including utrs). If equal to null, it is assumed that the annotations only contain CDS regions', default='null')
parser.add_argument("--rsq-type", type=str, help='Type of RNA-Seq alignment. Either RUM or GSNAP are currently supported', default='gsnap')
parser.add_argument("--cds-tag", type=str, help='The description tag as it appears in ANNOT_FILE for coding regions.', default='CDS')
parser.add_argument("--djob-num-nodes", type=int, help='A number greater or equal than 2 would allow the use of distribjob for parallelizing expensive tasks', default=40)
parser.add_argument("--dataset-id", type=str, help='RNA-Seq dataset id. Typically contains or makes reference to the organism and author names. For example "oocyst_booth"', default='oocyst_booth')
parser.add_argument("--sample-id", type=str, help='RNA-Seq sample id. Typically contains or makes reference to the condition, tissue, stage of development in which the reads were generated. For example "day3"', default='day0')
parser.add_argument("--rsq-orientation", type=str, help='Has a value of "N","F" or "R" if the RNA-Seq reads have orientation or not. For paired reads these values must be given in pairs, i.e. "NN", "FR" and so on', default="N")
line  = "";

try:
    args = parser.parse_args(sys.argv[1:])
    prop_dir =  args.PREFIX_PROPDIR+"_input";
    djob_dir = args.PREFIX_PROPDIR+"_DJob";
    log_file = open(args.PREFIX_PROPDIR+".error.log", "w")

    if os.path.exists(prop_dir) :
        shell_cmd = "rm -rf "+prop_dir
        ioutils.log_and_exec(shell_cmd, os.environ, log_file)

    if os.path.exists(djob_dir) :
        shell_cmd = "rm -rf "+djob_dir
        ioutils.log_and_exec(shell_cmd, os.environ, log_file)

    shell_cmd = "mkdir "+prop_dir
    ioutils.log_and_exec(shell_cmd, os.environ, log_file)
    shell_cmd = "mkdir "+djob_dir
    ioutils.log_and_exec(shell_cmd, os.environ, log_file)

    task_props = "RSqType="+args.rsq_type+"\nRSqOrientation="+args.rsq_orientation+"\nRSqWrapperIn=cat\nInputAnnotFormat=gff3\nInputContigFormat=fasta\ntranscriptTag="+args.transcript_tag+"\nCDSTag="+args.cds_tag+"\ngcClasses=100\ngeneModelType=ngscraig\nmodelUTRs=yes\ngenUTROnlyModel=yes\nnumDJobNodes="+str(args.djob_num_nodes)+"\nDJobNodeClass=SgeNode\nminPercJunctionAlignments=75\nspecies="+args.SPECIES+"\nRSqDataSetID="+args.dataset_id+"\nRSqSampleID="+args.sample_id+"\nRSqInputDir="+args.RSQ_INPUTDIR+"\nprefixOutputDir="+args.PREFIX_CRAIG_OUTPUT+"\nInputAnnotFile="+args.ANNOT_FILE+"\nInputContigFile=args.FASTA_FILE\n"

    controller_props = "masterdir="+prop_dir+"/master\ninputdir="+prop_dir+"\nnodedir="+djob_dir+"\nslotspernode=1\nsubtasksize=1\ntaskclass=DJob::DistribJobTasks::CRAIGTask\nnodeclass=DJob::DistribJob::SgeNode\n"
        
    taskprop_fd = open(prop_dir+"/task.prop", "w")
    taskprop_fd.write(task_props)
    taskprop_fd.close()
    controllerprop_fd = open(prop_dir+"/controller.prop", "w")
    controllerprop_fd.write(controller_props)
    controllerprop_fd.close()
except ioutils.CraigUndefined, msg:
    print msg
