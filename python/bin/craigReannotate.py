#!/usr/bin/env python
import struct, argparse, string, sys, numpy, os, warnings, subprocess, os.path, shutil, logging, re
sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils, ioutils

parser = argparse.ArgumentParser(description='Automated training/predicion wrapper tool for CRAIG v.1.1. It is a reannotation pipeline. It was mainly designed for eupathdb but it can be applied generally.\nWritten by Axel E. Bernal (abernal@seas.upenn.edu).');
parser.add_argument("SPECIES", type=str, help='Species name, given by the common suffix of files present in the models directory. These files specify the resources, filters and features used for constructing gene models and will be generated.', action='store')
parser.add_argument("ANNOT_FILE", type=str, help='Name of file with gene annotations', action='store')
parser.add_argument("FASTA_FILE", type=str, help='Name of file with fasta sequences or contigs', action='store')
parser.add_argument("-v", "--verbose", action='store_true', help='Turns on output of debugging information to stderr', dest='verbose')
parser.add_argument("--prefix-output-dir", type=str, help='Output directory for storing all output files.', default=".", dest='prefix_output_dir')
parser.add_argument("--model-utrs", action='store_true', help='Models and predicts UTRs\n', dest='model_utrs')
parser.add_argument("--pre-config", type=str, help='Name of file with preconfiguration information. See README for Details', default='', dest='preconfig_file')
parser.add_argument("--annot-fmt", type=str, help='Gene location format for ANNOT_FILE. One of locs, gff3, gtf or genbank. The format "locs" is an internal representation of genes. The format gff3 can also be used for handling the simpler gff format', default='gff3')
parser.add_argument("--contig-fmt", type=str, help='Gene location format for FASTA_FILE. One of fasta or genbank.', default='fasta', dest='contig_fmt')
parser.add_argument("--annot-tag", type=str, help='A unique identifier for the annotation\'s ANNOT_FILE source of origin. For example, sample/library/mapping tool/prediction tool', default='annot_tag')
parser.add_argument("--transcript-tag", type=str, help='The description tag as it appears in ANNOT_FILE for transcript regions(including utrs). If equal to null, it is assumed that the annotations only contain CDS regions', default='exon', dest='tr_tag')
parser.add_argument("--cds-tag", type=str, help='The description tag as it appears in ANNOT_FILE for coding regions.', default='CDS')
parser.add_argument("--gc-classes", type=str, help='GC isochore classes present in the genome. For example for human it would normally be 43,51,57', default='100')
parser.add_argument("--use-model-length-limits", action='store_true', help='Uses the genomic region length limits currently set in the input model files rather than adjusting them using the ranges found in the input existing annotation. Use this option when the input model files are up-to-date\n')
parser.add_argument("--utr-only-model", action='store_true', help='Produces a utr_only model in addition to the denovo one. GFF3 tags for each file will be gene:CRAIG_utr_only and gene:CRAIG_denovo -- for the gbrowse database\n', dest='utr_only_model')
parser.add_argument("--output-model-scores", action='store_true', help='Outputs gene structure scores for each structure component as computed by the model. The name of these output files will contain the keyword "gms"\n', dest='output_model_scores')
parser.add_argument("-f", "--force-train", action='store_true', help='Train model even if parameter files already exist', dest='force_train')
parser.add_argument("--model", type=str, help='The type of model to be learnt. Either craig for ab initio, ecraig for evidence integration or ngscraig for rna-seq driven gene finding', default='craig')
parser.add_argument("--closest-species", type=str, help='Closest species to be used if feature definitions for the model type to be learnt are not available for the current SPECIES', default="")
parser.add_argument("--min-perc-junc-aligns", type=int, help='Minimum percentage of genes with introns that have good junction alignments that are needed in order to use the junction information as the sole source for splicing');
parser.add_argument("--training-iterations", type=int, help='Number of training iterations needed for convergence', default=10)
parser.add_argument("--model-params", type=str, help='model parameters to use instead of performing training', default="")
parser.add_argument("--djob-num-nodes", type=int, help='A number greater or equal than 2 would allow the use of distribjob for parallelizing prediction', default=1)
parser.add_argument("--djob-input", type=str, help='distribjob input directory, i.e. where controller.prop will be located', default="")
parser.add_argument("--djob-node-class",type=str, help='distribjob node class, either SgeNode or LocalNode', default="LocalNode")
line  = "";

try:
# checking for sanity
    if 'CRAIG_HOME' not in os.environ:
        raise ioutils.CraigUndefined("CRAIG_HOME is not defined. Edit your bashrc file")

    my_env = os.environ
    craig_home = os.environ['CRAIG_HOME']
    my_env["PATH"] = craig_home+"/perl/bin:" + \
        craig_home+"/python/bin:" + \
        craig_home+"/bin:" + my_env.get("PATH", '')
    
    args = parser.parse_args(sys.argv[1:])

    # determining the model name
    model_name = args.SPECIES
    if args.model == 'ecraig':
        model_name += "-evid"
    elif args.model == 'ngscraig':
        model_name += '-rna'
    elif args.model != 'craig':
        raise ioutils.CraigUndefined("--model option undefined")

    model_output_dir = args.prefix_output_dir+".model"
    if not os.path.exists(model_output_dir):
        os.makedirs(model_output_dir)
    
    # opening log files
    log_file = open(model_output_dir+"/"+model_name+".error.log", "a")
    logging.basicConfig(filename=model_output_dir+"/"+model_name+".master.log", level=logging.DEBUG)

    #### PREPROCESSING ####
    preproc_output_dir = args.prefix_output_dir+".preproc" \
        if args.prefix_output_dir[0] == "/" \
        else os.path.dirname(os.getcwd())+"/"+args.prefix_output_dir+".preproc"

    shell_cmd = "craigPreprocess.py --config-file config"+(" --use-model-length-limits" if args.use_model_length_limits else "")+(" --pre-config-file "+args.preconfig_file if args.preconfig_file else "")+" --out-dir "+preproc_output_dir +(" --closest-species "+args.closest_species if args.closest_species != "" else "")+(" --annot-tag "+args.annot_tag if args.annot_tag != "" else "")+" --annot-fmt "+args.annot_fmt+" --contig-fmt "+args.contig_fmt+" --transcript-tag "+args.tr_tag+" --cds-tag "+args.cds_tag+" --gc-classes "+args.gc_classes+" --model "+args.model+" --num-permutations 1000 --block-length 20 "+("--model-utrs " if args.model_utrs else "")

    if args.djob_num_nodes > 1:
        shell_cmd += "--djob-num-nodes "+str(args.djob_num_nodes)+" --djob-node-class "+args.djob_node_class+" --djob-input "+args.djob_input+"/preprocess "

    shell_cmd += args.SPECIES+" "+args.ANNOT_FILE+" "+args.FASTA_FILE
    ioutils.log_and_exec(shell_cmd, my_env, log_file)

#### TRAINING and TESTING #####
    # Sanity check
    if not os.path.exists(preproc_output_dir+"/config"):
        raise ioutils.CraigUndefined("File "+preproc_output_dir+"/config needed for training but has not been created!.\nMost likely craigPreprocess.py crashed")

    # root prefix for all generated files
    model_prefix = craig_home+"/models/"+model_name
    shell_cmd = "grep -s PrefixGenomeFiles "+preproc_output_dir+"/config"
    prefix_genome = (subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split()[1]
    base_prefix_genome = os.path.basename(prefix_genome)
    parmodel_prefix = craig_home+"/models/"+base_prefix_genome

    shell_cmd = "grep -s Path "+preproc_output_dir+"/config"
    preproc_path = os.path.dirname((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split()[1])
    chrfasta_fn = preproc_path+"/"+model_name+".chr.fa"
    
    shell_cmd = "grep -s PrefixFiles "+preproc_output_dir+"/config"
    prefix_train = (subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split()[1]

    shell_cmd = "grep -s PrefixXValFiles "+preproc_output_dir+"/config"
    prefix_xvalfiles = (subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split()[1]

    shell_cmd = "cat "+prefix_xvalfiles+".locs"
    num_train_ids = len((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split("\n"))

    if num_train_ids < 300:
        logging.info('Insufficient number of tranining instances')
        exit(0)

    logging.info('Deciding whether to find splice signals or the ones provided by the junctions ')
    craigtrain_opts = ""
#train deciding whether to find splice signals or use the ones provided by the junctions
    if args.model == 'ngscraig':
        shell_cmd = "grep -s -v \"PAR_\" "+prefix_xvalfiles+".err | grep -s \"Added\" | perl -ne '$_ =~ /.+\s(\S+)$/; print \"$1\n\";' | sort -u | awk '{print \">\" $1}' > "+prefix_xvalfiles+".err.ids"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        shell_cmd = "cat "+prefix_xvalfiles+".err.ids"
        unsupp_ids = len((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split("\n"))
        shell_cmd = 'grep -s ">" '+prefix_xvalfiles+'.fa > '+prefix_xvalfiles+'.total.ids'
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        shell_cmd = "cat "+prefix_xvalfiles+".total.ids"
        total_ids = len((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split("\n"))
        logging.info('There are '+str(unsupp_ids)+' sequences with unsupported RNA-Seq evidence and '+str(total_ids)+' ids in total for training')
        craigtrain_opts += " --prefix-evidence="+prefix_train+" "
        
        if 100*unsupp_ids/total_ids < 100 - args.min_perc_junc_aligns:
            logging.info("RNASeq data set looks good, we will use junction information only")
            craigtrain_opts += "--DONOR-resource=RNAseq-signals --ACCEPTOR-resource=RNAseq-signals "
        else:
            logging.info("There are more than 25% of ids without junction support. We need to use all splice sites in the sequence for training\n")
#    elif args.model == 'craig' or args.model == 'ecraig':
        
# making sure all the needed files are present
    if not os.path.exists(model_prefix+".resources"):
        raise ioutils.CraigUndefined("File "+model_prefix+".resources needed for training but does not exist!\n")

    if not os.path.exists(model_prefix+".filters"):
        raise ioutils.CraigUndefined("File "+model_prefix+".filters needed for training but does not exist!\n")

    if not os.path.exists(model_prefix+".features"):
        raise ioutils.CraigUndefined("File "+model_prefix+".features needed for training but does not exist!\n")

    if not os.path.exists(model_prefix+".partial.top"):
        raise ioutils.CraigUndefined("File "+model_prefix+".partial.top needed for training but does not exist!\n")

    if not os.path.exists(model_prefix+".complete.top"):
        raise ioutils.CraigUndefined("File "+model_prefix+".complete.top needed for training but does not exist!\n")

    orig_mparams = parmodel_prefix+".params"
    mparams = orig_mparams
    new_mparams = orig_mparams
    utronly_mparams = parmodel_prefix+".utr_only.params"
    vcounter = 0
    while os.path.exists(new_mparams):
        mparams = new_mparams
        vcounter += 1
        new_mparams = orig_mparams+".v"+str(vcounter)

    if args.model_params != '':
        mparams = args.model_params

    if not os.path.exists(mparams) or args.force_train :
        shell_cmd = "cd "+model_output_dir+" && craigTrain "+craigtrain_opts+" --r=100 --model-utrs --loss-function=HAMMING --algorithm=ARROW -v --limit="+str(args.training_iterations)+" --avg-method=all "+preproc_output_dir+"/config params"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)    

#decide which iteration to pick for final parameter model
        shell_cmd = "grep -s Validation "+model_output_dir+"/"+model_name+".error.log"
        valid_lines = (subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split("\n")

        minimum = 100000000;
        ind_min = -1;
        for iteration in range(len(valid_lines) - 1):
            line = valid_lines[iteration]
#            print "\"", line, "\""
            tokens = line.split()

            this_min = float(tokens[3])+float(tokens[4])+40*(float(tokens[7])+float(tokens[8])) + 100*(float(tokens[10])+float(tokens[11]))
            if this_min < minimum :
                minimum = this_min
                ind_min = iteration
            
        if os.path.exists(mparams) :
            sys.stderr.write("File "+mparams+" exists! Creating a new version\n")
            mparams = new_mparams
            utronly_mparams = utronly_mparams+".v"+str(vcounter)
            sys.stderr.write("model parameters will be placed in "+mparams+"\nUser needs to purge old parameteres manually if needed\n")

        shell_cmd = "cp "+model_output_dir+"/params"+str(ind_min + 1)+" "+mparams
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

    output_file = model_output_dir+"/"+base_prefix_genome+".denovo"
    if args.djob_num_nodes > 1:
        add_opts = "predictUTRs="+("yes" if args.model_utrs else "no")+"\nformat=locs\n"
        ioutils.startCraigPredictDJob(args.djob_input+"/model",
                                      chrfasta_fn,
                                      mparams, prefix_genome+".chr",
                                      args.djob_num_nodes, 
                                      args.djob_node_class,
                                      add_opts, output_file+".chr.locs",
                                      my_env, log_file)
    else:
        shell_cmd = "craigPredict --prefix-evidence="+prefix_genome+ (" --predict-utrs " if args.model_utrs else " ")+" --format=locs "+mparams+" "+chrfasta_fn+" > "+output_file+".chr.locs"
        ioutils.log_and_exec(shell_cmd, my_env, log_file) 
        
    shell_cmd = "cat "+output_file+".chr.locs | locs2GTF.pl "+chrfasta_fn+" | gtf2gff3.pl | sed \'s/5UTR/UTR/g\' | sed \'s/3UTR/UTR/g\' | sed \'s/CRAIG/CRAIG_denovo/g\' > "+output_file+".chr.gff3"
    ioutils.log_and_exec(shell_cmd, my_env, log_file)

    if args.output_model_scores:
        #get denovo gene model scores
        if args.djob_num_nodes > 1:
            add_opts="haveUTRs=yes\n"
            ioutils.startScoreGenesDJob(args.djob_input+"/model",
                                        chrfasta_fn,
                                        output_file+".chr.locs",
                                        mparams,
                                        prefix_genome+".chr",
                                        args.djob_num_nodes,
                                        args.djob_node_class,
                                        add_opts,
                                        output_file+".chr.gms.locs",
                                        my_env, log_file)
        else:
            ioutils.produceSubContigAnnotation(chrfasta_fn,
                                               output_file+".chr.locs",
                                               output_file, True, True,
                                               my_env, log_file)
    
            shell_cmd = "score_genes "+("--have-utrs " if args.model_utrs else "")+"--prefix-evidence="+prefix_genome+".chr "+mparams+" "+output_file+".fa "+output_file+"locs > "+output_file+".gms.locs"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            shell_cmd = "moveLocsToOrigContigs.pl "+output_file+".subcontigs < "+output_file+".gms.locs > "+output_file+".chr.gms.locs"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

#edit params for utr-only
    if args.model == 'ngscraig' and args.model_utrs and args.utr_only_model:
        fdparams = open(parmodel_prefix+".utr_only.patch", "w")
        fdparams.write("__FILTER_SECTION__\nLazyFilter\tLog-RNASeq-CDSEdges\tEdgeAligner false SS TD\nlog(Filter)\tLog-RNASeq-CDSEdges-input\tFastaxFilter<EdgeAlignment,EdgeAlignment> RNAseq-edges Log-RNASeq-CDSEdges\n")
        fdparams.write("__FEATURE_SECTION__\nFeature\tLog-RNASeq-Exons-AlignType\tAlignmentType 1 Log-RNASeq-CDSEdges\nFeature\tLog-RNASeq-Exons-AlignWeight\tAlignmentWeight 1 Log-RNASeq-CDSEdges\nFeature\tLinearBin4Log-RNASeq-CDSWeights\tLinearFitBin<double> Log-RNASeq-Exons-AlignWeight 5 1.25 0 1\nFeature\tLinearBin4xLog-RNASeq-ExonWeight\tFeatureXFeature<UCHAR,double> Log-RNASeq-Exons-AlignType LinearBin4Log-RNASeq-CDSWeights -> CODING*\n")
        fdparams.write("__PARAM_SECTION__\nOverwrite\tStart-Background\t0\t0\t900256\nOverwrite\tStop-Background\t0\t0\t500842\nLinearBin4xLog-RNASeq-ExonWeight\t4\n1\t1\n0\t-0.655927\n0\t1\n5\t0\n-7.48136\n-1.9973\n0\n0\n-29.3232\n5\t0\n-3.48951\n-0.961739\n0\n0\n350.3071\n")
        fdparams.close()
        shell_cmd = "editCRAIGParamFile.pl "+parmodel_prefix+".utr_only.patch < "+mparams+" > "+utronly_mparams
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

        output_file = model_output_dir+"/"+base_prefix_genome+".utr_only"
        if args.djob_num_nodes > 1:
            add_opts = "predictUTRs=yes\nformat=locs\nSTART-resource=RNAseq-signals\nSTOP-resource=RNAseq-signals\n"
            ioutils.startCraigPredictDJob(args.djob_input+"/model",
                                          chrfasta_fn,
                                          utronly_mparams,
                                          prefix_genome+".chr.utr_only",
                                          args.djob_num_nodes,
                                          args.djob_node_class,
                                          add_opts, output_file+".chr.locs",
                                          my_env, log_file)
        else:
            shell_cmd = "craigPredict --prefix-evidence="+prefix_genome+".chr.utr_only --START-resource=RNAseq-signals --STOP-resource=RNAseq-signals --predict-utrs --format=locs "+utronly_mparams+" "+chrfasta_fn+" > "+output_file+".chr.locs"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

        shell_cmd = "cat "+output_file+".chr.locs | locs2GTF.pl "+chrfasta_fn+" | gtf2gff3.pl | sed \'s/5UTR/UTR/g\' | sed \'s/3UTR/UTR/g\' | sed \'s/CRAIG/CRAIG_utr_only/g\' > "+output_file+".chr.gff3"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

        if args.output_model_scores:
            if args.djob_num_nodes > 1:
                add_opts="haveUTRs=yes\n"
                ioutils.startScoreGenesDJob(args.djob_input+"/model",
                                            chrfasta_fn,
                                            output_file+".chr.locs",
                                            utronly_mparams,
                                            prefix_genome+".chr.utr_only",
                                            args.djob_num_nodes,
                                            args.djob_node_class,
                                            add_opts,
                                            output_file+".chr.gms.locs",
                                            my_env, log_file)
            else:
                ioutils.produceSubContigAnnotation(chrfasta_fn,
                                                   output_file+".chr.locs",
                                                   output_file, True, True,
                                                   my_env, log_file)

# get gene model scores for utr_only predictions
                shell_cmd = "score_genes --have-utrs --prefix-evidence="+prefix_genome+".chr.utr_only "+mparams+" "+output_file+".fa "+output_file+".locs > "+output_file+".gms.locs"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "moveLocsToOrigContigs.pl "+output_file+".subcontigs < "+output_file+".gms.locs > "+output_file+".chr.gms.locs"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)

    shell_cmd = "grep -v Broken "+model_output_dir+"/"+model_name+".error.log > "+model_output_dir+"/"+model_name+".error.log2; mv "+model_output_dir+"/"+model_name+".error.log2 "+model_output_dir+"/"+model_name+".error.log"
    ioutils.log_and_exec(shell_cmd, my_env, log_file)
except ioutils.CraigUndefined, msg:
    print msg
