#!/usr/bin/env python

import struct, argparse, string, sys, numpy, os, warnings, subprocess, os.path, shutil, logging, re, tempfile
sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils, ioutils
parser = argparse.ArgumentParser(description='Preprocessing script for CRAIG v.1.1, a discriminative gene prediction tool for eukarya.\nWritten by Axel E. Bernal (abernal@seas.upenn.edu).')

parser.add_argument("SPECIES", type=str, help='Species name, given by the common suffix of files present in the models directory. These files specify the resources, filters and features used for constructing gene models and will be generated.', action='store')
parser.add_argument("ANNOT_FILE", type=str, help='Name of file with gene annotations', action='store')
parser.add_argument("FASTA_FILE", type=str, help='Name of file with fasta sequences or contigs', action='store')
parser.add_argument("--config-file", type=str, help='Output configuration file to be used during training that has information about validation, training and where to find the file(s) with feature, filter and resource definitions. See README for details. If not specified, preprocessing will not setup the input sequences for training a model', default='')
parser.add_argument("-v", "--verbose", action='store_true', help='Turns on output of debugging information to stderr', dest='verbose')
parser.add_argument("--out-dir", type=str, help='Output directory for storing all output files.', default=".", dest='out_dir')
parser.add_argument("--pre-config-file", type=str, help='Name of file with preconfiguration information. See README for Details', default='')
parser.add_argument("--model", type=str, help='The type of model to be learnt. Either craig for ab initio, ecraig for evidence integration or ngscraig for rna-seq driven gene finding', default='ngscraig')
parser.add_argument("--annot-fmt", type=str, help='Gene location format for ANNOT_FILE. One of locs, gff3, gtf or genbank. The format "locs" is an internal representation of genes. The format gff3 can also be used for handling the simpler gff format', default='gff3')
parser.add_argument("--contig-fmt", type=str, help='Gene location format for FASTA_FILE. One of fasta or genbank.', default='fasta', dest='contig_fmt')
parser.add_argument("-l", "--file-list", action='store_true', help='Input files with gene annotations and fasta sequences are lists of filenames instead of flat files', dest='file_list')
parser.add_argument("-nx","--use-non-expressed", action='store_true', help='Use gene annotations that seem to be real sequence-wise, but have almost no RNA-Seq expression', dest='use_non_expressed_genes')
parser.add_argument("--annot-tag", type=str, help='A unique identifier for the annotation\'s ANNOT_FILE source of origin. For example, sample/library/mapping tool/prediction tool', default='annot_tag')
parser.add_argument("--transcript-tag", type=str, help='The description tag as it appears in ANNOT_FILE for transcript regions(including utrs). If equal to null, it is assumed that the annotations only contain CDS regions', default='exon', dest='tr_tag')
parser.add_argument("--cds-tag", type=str, help='The description tag as it appears in ANNOT_FILE for coding regions.', default='CDS')
parser.add_argument("--model-utrs", action='store_true', help='Models and predicts UTRs. Turned on when predicting NGS models\n', dest='model_utrs')
parser.add_argument("--use-default-length-limits", action='store_true', help='Uses the genomic region length limits currently set in the input model files rather than adjusting them using the ranges found in the input existing annotation. Use this option when the input model files are up-to-date')
parser.add_argument("-n", "--non-uniq-ids", action='store_true', help='Gene and transcript ids in ANNOT_FILE are not unique', dest='non_uniq_ids')
parser.add_argument("--fasta-wrapper", type=str, help='User-defined script executed before using FASTA_FILE if the input FASTA_FILE does not conform to the norm. Used in automated pipelines or settings', default='cat')
parser.add_argument("--annot-wrapper", type=str, help='User-defined script executed before using ANNOT_FILE if the input ANNOT_FILE does not conform to the norm. Used in automated pipelines or settings', default='cat')
parser.add_argument("--gc-classes", type=str, help='GC isochore classes present in the genome. For example for human it would normally be 43,51,57', default='100')
parser.add_argument("--min-coverage", type=int, help='Minium RNA-Seq coverage for a region transcription. It should be equal or greater than 1', default=2)
parser.add_argument("--smooth-window", type=int, help='Smoothing window for RNA-Seq  with non-uniform coverage for finding regions with active transcription.', default=0)
parser.add_argument("--num-permutations", type=int, help='Number of permutations to perform to detect change points in coverage', default=1000)
parser.add_argument("--closest-species", type=str, help='Closest species to be used if either an ab initio model or input feature definitions for the model type to be learnt are not available for the current SPECIES', default='generic')
parser.add_argument("--block-length", type=int, help='block length used in discarding permutations', default=20)
parser.add_argument("--djob-num-nodes", type=int, help='A number greater or equal than 2 would allow the use of distribjob for parallelizing expensive tasks', default=1)
parser.add_argument("--djob-input", type=str, help='distribjob input directory, i.e. where controller.prop will be located', default="")
parser.add_argument("--djob-node-class",type=str, help='distribjob node class, either SgeNode or LocalNode', default="LocalNode")
line  = ""

try:
# checking for sanity
    if 'CRAIG_HOME' not in os.environ:
        raise ioutils.CraigUndefined("CRAIG_HOME is not defined")
    my_env = os.environ
    craig_home = os.environ['CRAIG_HOME']
    my_env["PATH"] = craig_home+"/perl/bin:" + \
        craig_home+"/python/bin:" + \
        craig_home+"/bin:" + my_env.get("PATH", '')

    args = parser.parse_args(sys.argv[1:])

    if args.file_list:
        raise ioutils.CraigUndefined("--list-input option in development")

#    if os.path.exists(args.out_dir):
#        shutil.rmtree(args.out_dir)
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # determining the model name
    model_name = args.SPECIES
    if args.model == 'ecraig':
        model_name += "-evid"
    elif args.model == 'ngscraig':
        model_name += '-rna'
        args.model_utrs = True

    out_file = args.out_dir+"/"+model_name
    # opening log files
    log_file = open(out_file+".error.log", "w")
    logging.basicConfig(filename=out_file+".master.log",level=logging.DEBUG)

    chrfasta_fn = out_file+".chr.fa"
    chrsources_fn = out_file+".chr.sources"
    chrannot_fn = out_file+".chr.locs"
    chrlen_fn = out_file+".chr.len"
    chrlen2_fn = out_file+".chr.len2"
    chrss_fn = out_file+".chr.ssp"

    logging.info('Preparing fasta and annotation files for coverting to standard internal format')
    if not os.path.exists(out_file+".chr.prepfasta_is_done"):
        ioutils.processFasta(args.fasta_wrapper, args.contig_fmt, args.FASTA_FILE, out_file+".chr", my_env, log_file)
        ioutils.processAnnotation(args.annot_wrapper, args.annot_fmt, args.tr_tag, args.cds_tag, args.non_uniq_ids, args.ANNOT_FILE, out_file+".chr", my_env, log_file)
        ioutils.computePreprocAuxFiles(chrfasta_fn, chrannot_fn, chrlen_fn, chrlen2_fn, chrss_fn, my_env, log_file)
        if not args.model_utrs:
            shell_cmd = "stripUTRs.pl < "+chrannot_fn+" | filterRepeatedGenes.pl > "+chrannot_fn+".1"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            shell_cmd = "mv "+chrannot_fn+".1 "+chrannot_fn
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

        shell_cmd = "touch "+out_file+".chr.prepfasta_is_done"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

    # root prefix for all generated files
    model_prefix = craig_home+"/models/"+model_name
    naive_rsq_parmodel = out_file+".rsq_naive.params"
    if args.model == 'ngscraig':
        abi_model = craig_home+"/models/"+args.SPECIES+".params"
        if not os.path.exists(abi_model):
            abi_model = craig_home+"/models/"+args.closest_species+".params"

        if not os.path.exists(abi_model):
            raise ioutils.CraigUndefined("Cannot learn a ngscraig model type without an available ab initio model. Use option --closest-species correctly to provide a closely-related one")

        shell_cmd = "extractModelParams.pl -input "+abi_model+" -output "+out_file+".ab_initio.params < "+chrannot_fn
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        ioutils.makeNaiveRSqParamModel(naive_rsq_parmodel,
                                       out_file+".ab_initio.params",
                                       my_env, log_file)
        shell_cmd = "rm "+out_file+".ab_initio.params"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

    upd_length_limits = True

    if not os.path.exists(model_prefix+".features") :
        closest_model = args.closest_species;
        evid_sources = ""
        if args.model == 'ecraig':
            closest_model += "-evid"
            shell_cmd = "grep genepred "+args.pre_config_file+" | awk '{print $3 \"\_\" $4 \"\t\" $5}' > "+craig_home+"/models/"+model_name+".sources"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            evid_sources = " -sources "+craig_home+"/models/"+model_name+".sources"
        elif args.model == 'ngscraig':
            closest_model += '-rna'

        closest_model_path = craig_home+"/models/"+closest_model
        if not os.path.exists(closest_model_path+".template"):
            if not os.path.exists(closest_model_path+".features"):
                raise ioutils.CraigUndefined("Cannot build template for closely-related species "+closest_model)

            shell_cmd = "buildCRAIGModelTemplate.pl -output-template "+closest_model_path+".template -input-model "+closest_model_path
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

        shell_cmd = "buildCRAIGModelFilesFromTemplate.pl -input-template "+closest_model_path+".template -output-model "+model_prefix+" -extract-model-params "+chrannot_fn+evid_sources
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        upd_length_limits = False

    actual_model_name = model_name

    if not args.use_default_length_limits and upd_length_limits:
        shell_cmd = "buildCRAIGModelTemplate.pl -output-template "+out_file+".template -input-model "+model_prefix
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

        actual_model_name = out_file+".adjusted"
        shell_cmd = "buildCRAIGModelFilesFromTemplate.pl -input-template "+out_file+".template -output-model "+actual_model_name+" -extract-model-params "+chrannot_fn
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

    shell_cmd = "cat "+chrannot_fn
    num_ctglocs = len((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split("\n"))

    #generating either evidence or ab initio input files and features
    # pre_config_file contains information about sources of evidence
    # for rnaseq sources this is the header:
    # evid_type(rnaseq|genepred)  subtype(rum|gsnap|epdb*) dataset_id sample_id filename stranded(NN|N|F|R|FR|RF|FF|RR) wrapper_in

    prefix_genome = ""
    prefix_train = ""
    prefix_xvalfiles = ""

    # processing sources of evidence
    if len(args.pre_config_file):
        pre_config_fd = open(args.pre_config_file, "r")
        for line in pre_config_fd:
            line.strip()
            if line[0] == '#':
                continue
            tokens = line.split()
            logging.info('Processing evidence source '+tokens[0]+" "+tokens[4])
            rsqfn1 = out_file+"."+tokens[2]+"."+tokens[3]
            rsqfn2 = rsqfn1+".chr.junction"

            if tokens[0] == 'genepred':
                print "in development\n"
                if not len(prefix_genome):
                    prefix_genome = out_file+"."+args.annot_tag+".allev"
            elif tokens[0] == 'rnaseq':
                stranded = True if tokens[5].find('N') == -1 else False
                if len(prefix_genome) and prefix_genome != out_file+".allev":
                    raise ioutils.CraigUndefined("only one rnaseq source allowed")
                elif not len(prefix_genome):
                    prefix_genome = rsqfn1

                if tokens[1] == 'rum':
                    logging.info('Processing RNA-Seq reads for '+tokens[4])
                    if not os.path.exists(rsqfn1+".chr.rnaseq_is_done"):
                        shell_cmd = "sed 1,1d "+tokens[4]+"/junctions_all.rum | sort > "+out_file+".1sorted"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "sed 1,1d "+tokens[4]+"/junctions_high-quality.bed | perl -ne '@l = split(/\s+/, $_); $l[10] =~/([^\,]+)\,(\d+)/; $l[1] = $l[1] + $1 + 1; $l[2] = $l[2] - $2; print \"$l[0]\:$l[1]\-$l[2]\n\";' | sort > "+out_file+".2sorted"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "join "+out_file+".2sorted "+out_file+".1sorted | junctions2weightedLocs.py --orientation " + tokens[5] + " --dataset-id "+ tokens[2] + " --sample-id "+tokens[3]+" --cut-off 1 --format rum - " + rsqfn2 + ".orig.locs"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "join -v2 "+out_file+".2sorted "+out_file+".1sorted | junctions2weightedLocs.py --dataset-id " + tokens[2] + " --sample-id "+tokens[3] + " --cut-off 2 --format " + tokens[1] + " - - >> " + rsqfn2 + ".orig.locs"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "rm "+out_file+".1sorted "+out_file+".2sorted"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        logging.info('Processing coverage reads for '+tokens[4])

                        # dealing with coverage
                        if tokens[5] == 'N' or tokens[5] == 'NN':
                            shell_cmd = "sed 1,1d "+tokens[4]+"/RUM_Unique.cov | coverage2filter.py --len-file="+chrlen_fn+" --format=rum - > "+rsqfn1+".chr.cov"

                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        else:
                            cov_file = "RUM_Unique.plus.cov" if os.path.exists(tokens[4]+"/RUM_Unique.plus.cov") else "RUM_Unique_plus.cov"
                            shell_cmd = "sed 1,1d "+tokens[4]+"/"+cov_file+" | coverage2filter.py --len-file="+chrlen_fn+" --format=rum - | (perl -ne '$counter++; if($_ =~ /^>/ && $counter > 1) { print \"//\n\", $_; } else { print $_;}' && echo //)> "+rsqfn1+".chr.plus.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)

                            cov_file = "RUM_Unique.minus.cov" if os.path.exists(tokens[4]+"/RUM_Unique.minus.cov") else "RUM_Unique_minus.cov"
                            shell_cmd = "sed 1,1d "+tokens[4]+"/"+cov_file+" | coverage2filter.py --len-file="+chrlen_fn+" --format=rum - | (perl -ne '$counter++; if($_ =~ /^>/ && $counter > 1) { print \"//\n\", $_; } else { print $_;}' && echo //) > "+rsqfn1+".chr.minus.1.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            shell_cmd = "cat "+rsqfn1+".chr.minus.1.cov | moveScoreFilterToCompStrand.pl "+chrlen_fn +" > "+rsqfn1+".chr.minus.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            shell_cmd = "biosort.py "+rsqfn1+".chr.plus.cov > "+rsqfn1+".chr.plus.1.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            shell_cmd = "biosort.py "+rsqfn1+".chr.minus.cov > "+rsqfn1+".chr.minus.1.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            shell_cmd = "joinScoreFilters.pl "+rsqfn1+".chr.plus.1.cov "+rsqfn1+".chr.minus.1.cov > "+rsqfn1+".chr.cov"
                            if ioutils.log_and_exec(shell_cmd, my_env, log_file):
                                raise ioutils.CraigUndefined("joinScoreFilters.pl has found two different scores for one single position in input")

                            shell_cmd = "rm "+rsqfn1+".chr.plus.1.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            shell_cmd = "rm "+rsqfn1+".chr.minus.1.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        # add contigs without coverage
                        shell_cmd = "biodiff.py "+chrlen2_fn+" "+rsqfn1+".chr.cov -by ii >"+rsqfn1+".chr.cov.suppl ; cat "+rsqfn1+".chr.cov.suppl >> "+rsqfn1+".chr.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        shell_cmd = "touch "+rsqfn1+".chr.rnaseq_is_done"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                elif tokens[1] == 'gsnap':
                    if not os.path.exists(rsqfn1+".chr.rnaseq_is_done"):

                        #dealing with junctions
                        if not os.path.exists(tokens[4]+"/GSNAP.junction.locs"):
                            shell_cmd = "cat "+tokens[4]+"/GSNAP.aln | junctions2weightedLocs.py --dataset-id "+tokens[2]+" --orientation "+tokens[5]+" --format gsnap - - | filterRepeatedGenes.pl "+rsqfn2+".tmp.locs"
                        else:
                            shell_cmd = "cp "+tokens[4]+"/GSNAP.junction.locs "+rsqfn2+".tmp.locs"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        shell_cmd = "mv "+rsqfn2+".tmp.locs "+rsqfn2+".locs"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        #dealing with coverage
                        if not os.path.exists(tokens[4]+"/GSNAP.Unique.cov"):
                            shell_cmd = "cat "+tokens[4]+"/GSNAP.aln | gsreads2covfilter.py --orientation "+tokens[5]+" --len-file "+chrlen_fn+" > "+rsqfn1+".chr.cov"
                        else:
                            shell_cmd = "cp "+tokens[4]+"/GSNAP.Unique.cov "+rsqfn1+".chr.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        # add contigs without coverage
                        shell_cmd = "biodiff.py "+chrlen2_fn+" "+rsqfn1+".chr.cov -by ii >"+rsqfn1+".chr.cov.suppl ; cat "+rsqfn1+".chr.cov.suppl >> "+rsqfn1+".chr.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        shell_cmd = "touch "+rsqfn1+".chr.rnaseq_is_done"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                elif tokens[1] == 'bam':
                    if not os.path.exists(rsqfn1+".chr.rnaseq_is_done"):

                        #dealing with junctions
                        if not os.path.exists(rsqfn2+".locs"):
                            shell_cmd = "junctions2weightedLocs.py --with-bed " + rsqfn2 + ".bed --orientation " + tokens[5] + " --dataset-id "+tokens[2]+" --format bam " + tokens[4] + " - | filterRepeatedGenes.pl > " + rsqfn2 + ".tmp.locs"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            shell_cmd = "mv "+rsqfn2+".tmp.locs "+rsqfn2+".locs"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        #dealing with coverage
                        if not os.path.exists(rsqfn1+".chr.cov"):
                            shell_cmd = "coverage2filter.py --orientation " + tokens[5] + " --len-file "+chrlen_fn+" " + tokens[4] + " > "+rsqfn1+".tmp1.chr.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            # add contigs without coverage
                            shell_cmd = "biodiff.py "+chrlen2_fn+" "+rsqfn1+".tmp1.chr.cov -by ii >"+rsqfn1+".tmp2.chr.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            shell_cmd = "cat "+rsqfn1+".tmp1.chr.cov "+rsqfn1+".tmp2.chr.cov > "+rsqfn1+".chr.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)
                            shell_cmd = "rm -f " + rsqfn1 + ".tmp1.chr.cov " + rsqfn1 + ".tmp2.chr.cov"
                            ioutils.log_and_exec(shell_cmd, my_env, log_file)

                        shell_cmd = "touch "+rsqfn1+".chr.rnaseq_is_done"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                else:
                    raise ioutils.CraigUndefined("format "+tokens[1]+" not recognized")

                # dealing with junctions in non-stranded rnaseq sets
                logging.info('Processing junction reads strand information for '+tokens[4])

                if not os.path.exists(rsqfn2+".locs") or os.stat(rsqfn2+".locs").st_size == 0 :
                    shell_cmd = "cut -f1,2 "+rsqfn2+".orig.locs | filterNonCannonSS.pl -canon-pairs "+chrss_fn+" -truncated -ctg "+ chrfasta_fn+" -filter gene | grep -s \">\" | biointers.py "+rsqfn2+".orig.locs - > "+rsqfn2+".orig.fwd.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "cut -f1,2 "+rsqfn2+".orig.locs | reverseExons.pl | filterNonCannonSS.pl -canon-pairs "+chrss_fn+" -truncated -ctg "+chrfasta_fn+" -filter gene | grep -s \">\" | biointers.py "+rsqfn2+".orig.locs - | reverseExons.pl > "+rsqfn2+".orig.comp.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "cat "+rsqfn2+".orig.comp.locs "+rsqfn2+".orig.fwd.locs | biodiff.py "+rsqfn2+".orig.locs - -by ii > "+rsqfn2+".orig.unasigned.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "biointers.py "+rsqfn2+".orig.fwd.locs "+rsqfn2+".orig.comp.locs | biodiff.py "+rsqfn2+".orig.fwd.locs - > "+rsqfn2+".locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "biointers.py "+rsqfn2+".orig.comp.locs "+rsqfn2+".orig.fwd.locs | biodiff.py "+rsqfn2+".orig.comp.locs - >> "+rsqfn2+".locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

                # compute junction.signal and .xsignal for rsqfn1+".chr" prefix
                if not os.path.exists(rsqfn1+".chr.xsig_is_done"):
                    logging.info('Computing junction signal filter for '+tokens[4])
                    shell_cmd = "extractWeightedSignalFilter.pl -ctg "+chrfasta_fn+" -truncated -locs "+rsqfn2
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    if args.djob_num_nodes > 1:
                        add_opts="minCoverage=2\nblockLen=20\nstranded="+("yes" if stranded else "no")+"\nnumPermutations="+str(args.num_permutations)+"\n"
                        ioutils.startXSignalsDJob(args.djob_input,
                                                  chrfasta_fn,
                                                  rsqfn1+".chr",
                                                  args.djob_num_nodes,
                                                  args.djob_node_class,
                                                  add_opts,
                                                  my_env, log_file)
                    else:
                        shell_cmd = "get_covxsignals "+("--stranded" if stranded else "")+" --prefix-evidence="+rsqfn1+".chr --num-permutations="+str(args.num_permutations)+" --block-len=20 --chunk-size=30000 "+chrfasta_fn
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                    shell_cmd = "touch "+rsqfn1+".chr.xsig_is_done"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

                #compute junction.signal and .xsignal for rsqfn1+".chr.utr_only
                #prefix
                if not os.path.exists(rsqfn1+".chr.utr_only.xsig_is_done") :
                    shell_cmd = "stripUTRs.pl < "+chrannot_fn+" | filterRepeatedGenes.pl > "+out_file+".chr.nutr.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "filterOlapGenes.pl -o "+out_file+".chr.nutr.locs < "+rsqfn2+".locs > "+rsqfn1+".chr.utr_only.junction.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "extractSignalFilter.pl --truncated -ctg "+chrfasta_fn+" < "+rsqfn1+".chr.utr_only.junction.locs > "+rsqfn1+".chr.utr_only.junction.1.signal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "cat "+out_file+".chr.nutr.locs >> "+rsqfn1+".chr.utr_only.junction.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "extractSignalFilter.pl --sig-weight=100 -ctg "+chrfasta_fn+" < "+out_file+".chr.nutr.locs > "+rsqfn1+".chr.utr_only.junction.2.signal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "biosort.py "+rsqfn1+".chr.utr_only.junction.1.signal > "+rsqfn1+".chr.utr_only.junction.3.signal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "biosort.py "+rsqfn1+".chr.utr_only.junction.2.signal > "+rsqfn1+".chr.utr_only.junction.4.signal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "joinSignalFilters.pl "+rsqfn1+".chr.utr_only.junction.4.signal "+rsqfn1+".chr.utr_only.junction.3.signal > "+rsqfn1+".chr.utr_only.junction.signal"
                    if ioutils.log_and_exec(shell_cmd, my_env, log_file):
                        raise ioutils.CraigUndefined("Execution of joinSignalFilters.pl has failed. Either two different signals for one single position were found, or some unknown mistake")
                    # link utr_only model files
                    shell_cmd = "ln -fs "+rsqfn1+".chr.cov "+rsqfn1+".chr.utr_only.cov"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "ln -fs "+rsqfn1+".chr.xsignal "+rsqfn1+".chr.utr_only.xsignal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "ln -fs "+rsqfn1+".chr.xsigscores "+rsqfn1+".chr.utr_only.xsigscores"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

                    shell_cmd = "touch "+rsqfn1+".chr.utr_only.xsig_is_done"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

                # computing tentative subcontig regions and filters to
                # extract single-gene sequences for training
                logging.info('Computing subcontig regions and filters')

                if not os.path.exists(rsqfn1+".subcontig_is_done") :
                    #finding intergenic regions
                    if args.djob_num_nodes > 1:
                        add_opts = "predictUTRs=no\nformat=locs\nDONOR-resource=RNAseq-signals\nACCEPTOR-resource=RNAseq-signals\n"
                        ioutils.startCraigPredictDJob(args.djob_input,
                                                      chrfasta_fn,
                                                      naive_rsq_parmodel,
                                                      rsqfn1+".chr",
                                                      args.djob_num_nodes,
                                                      args.djob_node_class,
                                                      add_opts,
                                                      rsqfn1+".chr.locs",
                                                      my_env, log_file)
                    else:
                        shell_cmd = "craigPredict --prefix-evidence="+rsqfn1+".chr --format=locs --DONOR-resource=RNAseq-signals --ACCEPTOR-resource=RNAseq-signals "+naive_rsq_parmodel+" "+chrfasta_fn+" > "+rsqfn1+".chr.locs"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

                    # using the annotation above, find the corresponding subcontigs
                    ioutils.produceSubContigAnnotation(chrfasta_fn,
                                                       rsqfn1+".chr.locs",
                                                       rsqfn1, False, False,
                                                       my_env, log_file)

                    # compute subcontig filters
                    ioutils.computeSubContigFilters(rsqfn1+".chr", rsqfn1,
                                                    args.djob_input,
                                                    args.djob_num_nodes,
                                                    args.djob_node_class,
                                                    args.num_permutations,
                                                    stranded,
                                                    my_env, log_file)

                    shell_cmd = "touch "+rsqfn1+".subcontig_is_done"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

                # compute utr_only filters for the subcontigs above
                if not os.path.exists(rsqfn1+".utr_only.xsig_is_done"):
                    shell_cmd = "stripUTRs.pl < "+rsqfn1+".locs | filterRepeatedGenes.pl > "+rsqfn1+".nutr.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "filterOlapGenes.pl -o "+rsqfn1+".nutr.locs < "+rsqfn1+".junction.locs > "+rsqfn1+".utr_only.junction.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "extractSignalFilter.pl --truncated -ctg "+rsqfn1+".fa < "+rsqfn1+".utr_only.junction.locs > "+rsqfn1+".utr_only.junction.1.signal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "cat "+rsqfn1+".nutr.locs >> "+rsqfn1+".utr_only.junction.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "extractSignalFilter.pl --sig-weight=100 -ctg "+rsqfn1+".fa < "+rsqfn1+".nutr.locs > "+rsqfn1+".utr_only.junction.2.signal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "biosort.py "+rsqfn1+".utr_only.junction.1.signal > "+rsqfn1+".utr_only.junction.3.signal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "biosort.py "+rsqfn1+".utr_only.junction.2.signal > "+rsqfn1+".utr_only.junction.4.signal"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "joinSignalFilters.pl "+rsqfn1+".utr_only.junction.4.signal "+rsqfn1+".utr_only.junction.3.signal > "+rsqfn1+".utr_only.junction.signal"
                    if ioutils.log_and_exec(shell_cmd, my_env, log_file):
                        raise ioutils.CraigUndefined("joinSignalFilters.pl has found two different signals for one single position in input")

                    shell_cmd = "touch "+rsqfn1+".utr_only.xsig_is_done"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

            else:
                raise ioutils.CraigUndefined("Unrecognized source type")

    if not len(prefix_genome) and args.model == 'craig':
        # compute ab initio gene split of the annotation
        prefix_genome = out_file+"."+args.annot_tag+".abi"
        if not os.path.exists(prefix_genome+".final.subset_is_done"):
            shell_cmd = "cat "+chrannot_fn+" | extractIGenics.pl -ctg "+chrfasta_fn+" > "+prefix_genome+".igenics"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            ioutils.produceSubsetAnnotation(chrfasta_fn,
                                            prefix_genome+".igenics",
                                            chrannot_fn,
                                            prefix_genome,
                                            my_env, log_file)

            shell_cmd = "genes2tags "+("--have-utrs " if args.model_utrs else "")+"--species="+actual_model_name+" "+prefix_genome+".fa "+prefix_genome+".locs 2> "+prefix_genome+".err > "+prefix_genome+".tag"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            shell_cmd = "grep -s -v \"START_\" < "+prefix_genome+".err | grep -s -v \"STOP_\" | grep -s \"Added\" | perl -ne '$_ =~ /.+\s(\S+)$/; print \"$1\n\";' | sort -u | awk '{print \">\" $1}' > "+prefix_genome+".defect.ids"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            shell_cmd = "biodiff.py "+prefix_genome+".locs "+prefix_genome+".defect.ids -by ci > "+prefix_genome+".final.locs"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            shell_cmd = "biodiff.py "+prefix_genome+".fa "+prefix_genome+".defect.ids -by ii > "+prefix_genome+".final.fa"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            shell_cmd = "biodiff.py "+prefix_genome+".tag "+prefix_genome+".defect.ids -by ii > "+prefix_genome+".final.tag"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            shell_cmd = "touch "+prefix_genome+".final.subset_is_done"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

        prefix_train = prefix_genome+".final"

    elif not len(prefix_genome):
        raise ioutils.CraigUndefined("empty pre-configuration file or --model option undefined")

    # obtaining a training set
    if args.model == 'ngscraig':
        if not os.path.exists(prefix_genome+".trainingset_is_done"):
            if args.djob_num_nodes > 1:
                add_opts = "predictUTRs=no\nformat=locs\nDONOR-resource=RNAseq-signals\nACCEPTOR-resource=RNAseq-signals\n"
                ioutils.startCraigPredictDJob(args.djob_input,
                                              prefix_genome+".fa",
                                              naive_rsq_parmodel,
                                              prefix_genome,
                                              args.djob_num_nodes,
                                              args.djob_node_class,
                                              add_opts,
                                              prefix_genome+".1.locs",
                                              my_env, log_file)
            else:
                shell_cmd = "craigPredict --prefix-evidence="+prefix_genome+" --format=locs --DONOR-resource=RNAseq-signals --ACCEPTOR-resource=RNAseq-signals "+naive_rsq_parmodel+" "+prefix_genome+".fa > "+prefix_genome+".1.locs"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)

            shell_cmd = "cat "+prefix_genome+".nutr.locs | grep -s -v \"<\" | grep -E -s -v \"\S>\" | filterOlapGenes.pl -v -o "+prefix_genome+".1.locs > "+prefix_genome+".2.locs"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            shell_cmd = "remove_suspect_transcripts --prefix-evidence="+prefix_genome+" "+prefix_genome+".fa "+prefix_genome+".nutr.locs "+prefix_genome+".2.locs | biointers.py -by ii "+prefix_genome+".nutr.locs - > "+prefix_genome+".3.locs"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

            if not args.use_non_expressed_genes:
                shell_cmd = "cat "+prefix_genome+".3.locs  | filterOlapGenes.pl -s --min-weight 2 -v -o "+prefix_genome+".junction.locs > "+prefix_genome+".4.locs"
            else:
                shell_cmd = "mv "+prefix_genome+".3.locs "+prefix_genome+".4.locs"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

            shell_cmd = "touch "+prefix_genome+".trainingset_is_done"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

        shell_cmd ="removeGenes.sh "+prefix_genome+".nutr.locs "+prefix_genome+".4.locs > "+prefix_genome+".naivPrCorr.noAS.locs"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        shell_cmd = "cat "+prefix_genome+".naivPrCorr.noAS.locs"
        num_train_ids = len((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = my_env, stderr = log_file).communicate()[0]).split("\n"))
        logging.info('There are '+str(num_train_ids)+' sequences for training')

        if num_train_ids > 0.2*num_ctglocs or num_train_ids > 4000:
            shell_cmd = "biointers.py "+prefix_genome+".fa "+prefix_genome+".naivPrCorr.noAS.locs -by ic > "+prefix_genome+".naivPrCorr.noAS.fa"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)
            prefix_train = prefix_genome
        else:
            prefix_output = prefix_genome+".naivPrCorr.noAS"
            if not os.path.exists(prefix_output+".trainingset_is_done"):
                shell_cmd = "biointers.py -by cc "+prefix_genome+".1.locs "+prefix_genome+".4.locs | extractIGenics.pl -ctg "+prefix_genome+".fa > "+prefix_output+".1.igenics"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "biointers.py -by cc "+rsqfn1+".nutr.locs "+prefix_genome+".4.locs | extractIGenics.pl -ctg "+prefix_genome+".fa | filterOlapGenes.pl -o "+prefix_output+".1.igenics > "+prefix_output+".2.igenics"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "cat "+prefix_output+".1.igenics "+prefix_output+".2.igenics > "+prefix_output+".igenics"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)

                ioutils.produceSubsetAnnotation(prefix_genome+".fa",
                                                prefix_output+".igenics",
                                                prefix_genome+".4.locs",
                                                prefix_output,
                                                my_env, log_file)

                ioutils.computeSubContigFilters(prefix_genome,
                                                prefix_output,
                                                args.djob_input,
                                                args.djob_num_nodes,
                                                args.djob_node_class,
                                                args.num_permutations,
                                                stranded,
                                                my_env, log_file)

                shell_cmd = "touch "+prefix_output+".trainingset_is_done"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)

            prefix_train = prefix_genome+".naivPrCorr.noAS"

        if not os.path.exists(prefix_genome+".training_utrannot_is_done"):
            ioutils.annotateUTRs(prefix_train,
                                 prefix_genome+".naivPrCorr.noAS",
                                 actual_model_name, stranded, my_env, log_file)
            shell_cmd = "touch "+prefix_genome+".training_utrannot_is_done"
            ioutils.log_and_exec(shell_cmd, my_env, log_file)

# for now we will proceed without alternative splicing filtering
#        shell_cmd = "nohup report_rnaseq_asevents --prefix-evidence="+prefix_fil#es+" "+prefix_genome+".fa > "+prefix_genome+".ASreport"
#        ioutils.log_and_exec(shell_cmd, my_env, log_file)
#        shell_cmd = "grep -s \"is retained\" "+prefix_genome+".ASreport | cut -f1 -d \" \" > "+prefix_genome+".6.locs"
#        ioutils.log_and_exec(shell_cmd, my_env, log_file)
#        shell_cmd = "grep -s \"exon skip\" "+prefix_genome+".ASreport | cut -f1 -d \" \" | cat - "+prefix_genome+".6.locs | sort -u | perl -ne 'print \">$_\";' | sort -u | biodiff.py "+prefix_genome+".naivPrCorr.locs - -by ci > "+prefix_genome+".naivPrCorr.noAS.locs"
#        ioutils.log_and_exec(shell_cmd, my_env, log_file)

        prefix_xvalfiles = prefix_genome+".naivPrCorr.noAS.utr"

    elif args.model == 'ecraig':
        prefix_xvalfiles = prefix_genome
    elif args.model == 'craig':
        prefix_xvalfiles = prefix_genome+".final"
    else:
        raise ioutils.CraigUndefined("--model option undefined")

    if len(args.config_file):
        ioutils.createSetup4Training(prefix_train, prefix_genome,
                                     prefix_xvalfiles,
                                     actual_model_name, args.gc_classes,
                                     args.out_dir, args.config_file,
                                     my_env, log_file)

except ioutils.CraigUndefined, msg:
    print msg
