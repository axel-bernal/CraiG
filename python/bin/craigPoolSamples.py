#!/usr/bin/env python

import struct, argparse, string, sys, numpy, os, warnings, subprocess, os.path, shutil, logging, re
sys.path.append(os.environ['CRAIG_HOME']+'/python/lib/')

import bioutils, ioutils
parser = argparse.ArgumentParser(description='Preprocessing script for CRAIG v.1.1, a discriminative gene prediction tool for eukarya.\nWritten by Axel E. Bernal (abernal@seas.upenn.edu).')

parser.add_argument("SPECIES", type=str, help='Species name, given by the common suffix of files present in the models directory. These files specify the resources, filters and features used for constructing gene models and will be generated.', action='store')
parser.add_argument("ANNOT_FILE", type=str, help='Name of file with gene annotations', action='store')
parser.add_argument("FASTA_FILE", type=str, help='Name of file with fasta sequences or contigs', action='store')
parser.add_argument("CONF_FILE", type=str, help='Output configuration file to be used during training that has information about validation, training and where to find the file(s) with feature, filter and resource definitions. See README for details')
parser.add_argument("-v", "--verbose", action='store_true', help='Turns on output of debugging information to stderr', dest='verbose')
parser.add_argument("--out-dir", type=str, help='Output directory for storing all output files.', default=".", dest='out_dir')
parser.add_argument("--pre-config", type=str, help='Name of file with preconfiguration information. See README for Details', default='', dest='preconfig_file')
parser.add_argument("--model", type=str, help='The type of model to be learnt. Either craig for ab initio, ecraig for evidence integration or ngscraig for rna-seq driven gene finding', default='craig')
parser.add_argument("--ab-initio-model", type=str, help='ab initio model needed if the model to compute is either ecraig or ngscraig', default='')
parser.add_argument("--annot-fmt", type=str, help='Gene location format for ANNOT_FILE. One of locs, gff3, gtf or genbank. The format "locs" is an internal representation of genes. The format gff3 can also be used for handling the simpler gff format', default='gff3')
parser.add_argument("--contig-fmt", type=str, help='Gene location format for FASTA_FILE. One of fasta or genbank.', default='fasta', dest='contig_fmt')
parser.add_argument("-l", "--file-list", action='store_true', help='Input files with gene annotations and fasta sequences are lists of filenames instead of flat files', dest='file_list')
parser.add_argument("-nx","--use-non-expressed", action='store_true', help='Use gene annotations that seem to be real sequence-wise, but have almost no RNA-Seq expression', dest='use_non_expressed_genes')
parser.add_argument("--transcript-tag", type=str, help='The description tag as it appears in ANNOT_FILE for transcript regions(including utrs). If equal to null, it is assumed that the annotations only contain CDS regions', default='null', dest='tr_tag')
parser.add_argument("--cds-tag", type=str, help='The description tag as it appears in ANNOT_FILE for coding regions.', default='CDS')
parser.add_argument("-n", "--non-uniq-ids", action='store_true', help='Gene and transcript ids in ANNOT_FILE are not unique', dest='non_uniq_ids')
parser.add_argument("--fasta-wrapper", type=str, help='User-defined script executed before using FASTA_FILE if the input FASTA_FILE does not conform to the norm. Used in automated pipelines or settings', default='cat')
parser.add_argument("--annot-wrapper", type=str, help='User-defined script executed before using ANNOT_FILE if the input ANNOT_FILE does not conform to the norm. Used in automated pipelines or settings', default='cat')
parser.add_argument("--gc-classes", type=str, help='GC isochore classes present in the genome. For example for human it would normally be 43,51,57', default='100')
parser.add_argument("--min-coverage", type=int, help='Minium RNA-Seq coverage for a region transcription. It should be equal or greater than 1', default=2)
parser.add_argument("--coverage-density", type=int, help='The coverage density factor. A number between 1 (very sparse) to 9 (very dense)', default=5)
parser.add_argument("--smooth-window", type=int, help='Smoothing window for RNA-Seq  with non-uniform coverage for finding regions with active transcription.', default=0)
parser.add_argument("--num-permutations", type=int, help='Number of permutations to perform to detect change points in coverage', default=1000)
parser.add_argument("--block-length", type=int, help='block length used in discarding permutations', default=20)
parser.add_argument("--djob-num-nodes", type=int, help='A number greater or equal than 2 would allow the use of distribjob for parallelizing expensive tasks', default=1)
parser.add_argument("--djob-input", type=str, help='distribjob input directory, i.e. where controller.prop will be located', default="")
parser.add_argument("--djob-node-class",type=str, help='distribjob node class, either SgeNode or LocalNode', default="LocalNode")
line  = "";

class CraigUndefined(Exception):
    pass

try:
# checking for sanity
    if 'CRAIG_HOME' not in os.environ:
        raise CraigUndefined("CRAIG_HOME is not defined")
    my_env = os.environ
    my_env["PATH"] = os.environ['CRAIG_HOME']+"/perl/bin:" + \
        os.environ['CRAIG_HOME']+"/python/bin:" + \
        os.environ['CRAIG_HOME']+"/bin:" + my_env.get("PATH", '')
    
    args = parser.parse_args(sys.argv[1:])

    if args.file_list:
        raise CraigUndefined("--list-input option in development")

    if os.path.exists(args.out_dir):
        shutil.rmtree(args.out_dir)

    os.makedirs(args.out_dir)
    
    # determining the model name
    model_name = args.SPECIES
    if args.model == 'ecraig':
        model_name += "-evid"
    elif args.model == 'ngscraig':
        model_name += '-rna'

    outFile = args.out_dir+"/"+model_name
    # opening log files
    log_file = open(outFile+".error.log", "w")
    logging.basicConfig(filename=outFile+".master.log",level=logging.DEBUG)

    # root prefix for all generated files
    model_prefix = os.environ['CRAIG_HOME']+"/models/"+model_name
    ab_model = args.ab_initio_model
    if ab_model == '':
        ab_model = os.environ['CRAIG_HOME']+"/models/"+args.SPECIES+".params"
    if args.model == 'ngscraig' and not os.path.exists(ab_model):
        raise CraigUndefined("ab initio model "+ab_model+" must be available for RNA-Seq")
    
    logging.info('Transforming annotation and fasta to standard internal format')
    ioutils.processFasta(args.fasta_wrapper, args.contig_fmt, args.FASTA_FILE, outFile+".chr", my_env, log_file)
    ioutils.processAnnotation(args.annot_wrapper, args.annot_fmt, args.tr_tag, args.cds_tag, args.non_uniq_ids, args.ANNOT_FILE, outFile+".chr", my_env, log_file)

    chrfasta_fn = outFile+".chr.fa"
    chrsources_fn = outFile+".chr.sources"
    chrannot_fn = outFile+".chr.locs"
    chrlen_fn = outFile+".chr.len"
    chrlen2_fn = outFile+".chr.len2"
    chrss_fn = outFile+".chr.ssp"
    logging.info('Computing auxiliary files')
    ioutils.computePreprocAuxFiles(chrfasta_fn, chrannot_fn, chrlen_fn, chrlen2_fn, chrss_fn, my_env, log_file)

    #generating either evidence or ab initio input files and features
    # preconfig file contains information about sources of evidence
    # for rnaseq sources this is the header:
    # evid_type(rnaseq|gpred)  subtype(rum|gsnap|epdb*) dataset_id sample_id filename stranded(yes|no) wrapper_in
    
    prefix_files = ""
    prefix_trainfiles = ""

    # processing sources of evidence
    if len(args.preconfig_file):
        preconfig_fd = open(args.preconfig_file, "r")
        for line in preconfig_fd:
            line.strip()
            if line[0] == '#':
                continue
            tokens = line.split()
            logging.info('Processing evidence source '+tokens[0]+" "+tokens[4])
            rsqfn1 = outFile+"."+tokens[2]+"."+tokens[3]
            rsqfn2 = rsqfn1+".chr.junction";

            if tokens[0] == 'gpred':
                print "in development\n"
                if not len(prefix_files):
                    prefix_files = outFile+".allev"
            elif tokens[0] == 'rnaseq':
                if len(prefix_files) and prefix_files != outFile+".allev":
                    raise CraigUndefined("only one rnaseq source allowed")
                elif not len(prefix_files):
                    prefix_files = rsqfn1

                if tokens[1] == 'rum':
                    logging.info('Processing junction reads for '+tokens[4])
                    # dealing with junctions
                    shell_cmd = "sed 1,1d "+tokens[4]+"/junctions_all.rum | sort > "+outFile+".1sorted" 
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "sed 1,1d "+tokens[4]+"/junctions_high-quality.bed | perl -ne '@l = split(/\s+/, $_); $l[10] =~/([^\,]+)\,(\d+)/; $l[1] = $l[1] + $1 + 1; $l[2] = $l[2] - $2; print \"$l[0]\:$l[1]\-$l[2]\n\";' | sort > "+outFile+".2sorted"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "join "+outFile+".2sorted "+outFile+".1sorted | junctions2weightedLocs.py "+("--stranded" if tokens[5] == 'yes' else "")+" --dataset-id "+ tokens[2] + " --sample-id "+tokens[3]+" --cut-off 1 --format rum > "+rsqfn2+".orig.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "join -v2 "+outFile+".2sorted "+outFile+".1sorted | junctions2weightedLocs.py --dataset-id "+tokens[2] +" --sample-id "+tokens[3]+" --cut-off 2 --format "+tokens[1]+" >> "+rsqfn2+".orig.locs"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    shell_cmd = "rm "+outFile+".1sorted "+outFile+".2sorted"                   
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

                    logging.info('Processing coverage reads for '+tokens[4])                    
                     # dealing with coverage
                    if tokens[5] == 'no':
                        shell_cmd = "sed 1,1d "+tokens[4]+"/RUM_Unique.cov | coverage2filter.py --len-file="+chrlen_fn+" --format=rum > "+rsqfn1+".chr.cov"
 # add contigs without coverage
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    else:
                        shell_cmd = "sed 1,1d "+tokens[4]+"/RUM_Unique.plus.cov | coverage2filter.py --len-file="+chrlen_fn+" --format=rum | (perl -ne '$counter++; if($_ =~ /^>/\
 && $counter > 1) { print \"//\n\", $_; } else { print $_;}' && echo //)> "+rsqfn1+".chr.plus.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "sed 1,1d "+tokens[4]+"/RUM_Unique.minus.cov | coverage2filter.py --len-file="+chrlen_fn+" --format=rum | (perl -ne '$counter++; if($_ =~ /^>/ && $counter > 1) { print \"//\n\", $_; } else { print $_;}' && echo //) | moveScoreFilterToCompStrand.pl "+chrfasta_fn +" > "+rsqfn1+".chr.minus.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "biosort.py "+rsqfn1+".chr.plus.cov > "+rsqfn1+".chr.plus.1.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "biosort.py "+rsqfn1+".chr.minus.cov > "+rsqfn1+".chr.minus.1.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "joinScoreFilters.pl "+rsqfn1+".chr.plus.1.cov "+rsqfn1+".chr.minus.1.cov > "+rsqfn1+".chr.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        if not ioutils.checkLog4Completion("Done", outFile+".error.log", my_env, log_file):
                            raise CraigUndefined("joinScoreFilters.pl has found two different scores for one single position in input")
                        shell_cmd = "rm "+rsqfn1+".chr.plus.1.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)
                        shell_cmd = "rm "+rsqfn1+".chr.minus.1.cov"
                        ioutils.log_and_exec(shell_cmd, my_env, log_file)

 # add contigs without coverage
                    shell_cmd = "biodiff.py "+chrlen2_fn+" "+rsqfn1+".chr.cov -by ii >"+rsqfn1+".chr.cov.suppl ; cat "+rsqfn1+".chr.cov.suppl >> "+rsqfn1+".chr.cov"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)       
                        
#                    print 'hola'
                elif tokens[1] == 'gsnap':
                    #dealing with junctions
                    shell_cmd = "cat "+tokens[4]+"/\*concordant_uniq | junctions2weightedLocs.py "+("--stranded" if tokens[5] == 'yes' else "")+" --dataset-id "+tokens[2]+" --format gsnap | filterRepeatedGenes.pl perl -ne 'if($_ =~ /^>\S+\.(\d)\s+\S+\s+(\S+)/) { if($1 eq '0' && $2 eq \"1;1\") {;} else {print $_;}}' > "+rsqfn1+".chr.cov"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    #dealing with coverage
                    shell_cmd = "cat "+tokens[4]+"/\*concordat_uniq | gsreads2covfilter.py "+("--stranded" if tokens[5] == 'yes' else "")+" --len-file "+chrlen_fn+" > "+rsqfn1+".chr.cov"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                else:
                    shell_cmd = "sed 1,1d "+tokens[4]+" | junctions2weightedLocs.py --dataset-id "+ tokens[2] + "--sample-id "+tokens[3]+" --cut-off 1 --format "+tokens[1]
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

                # dealing with junctions in non-stranded rnaseq sets
                logging.info('Processing junction reads strand information for '+tokens[4])
#                if tokens[5] == 'no':
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
#                print 'hola'
#                else:
#                    shell_cmd = "mv "+rsqfn2+".orig.locs "+rsqfn2+".locs"
#                    ioutils.log_and_exec(shell_cmd, my_env, log_file)

# computing junction signal info for chromosomal level
                logging.info('Computing junction signal filter for '+tokens[4])
                shell_cmd = "extractWeightedSignalFilter.pl -ctg "+chrfasta_fn+" -truncated -locs "+rsqfn2
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "perl -ne '$_ =~ /^(\S+)\s+(\S+)\s*$/; print \">$2\t$1\t3\n\- x $1\n\/\/\n\- x $1\n\";' < "+chrlen_fn+" > "+rsqfn1+".chr.xsignal"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "perl -ne '$_ =~ /^(\S+)\s+(\S+)\s*$/; print \">$2\t$1\t3\n\/\/\n\";' < "+chrlen_fn+" > "+rsqfn1+".chr.xsigscores"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)

# moving information to subcontigs
                shell_cmd = "find_nocov_regions --coverage-density="+str(args.coverage_density)+" --max-genic=100000 --smooth-window="+str(args.smooth_window)+(" --stranded" if tokens[5] == 'yes' else " ")+ " --min-coverage="+str(args.min_coverage)+" --prefix-evidence="+rsqfn1+".chr "+chrfasta_fn+" > "+rsqfn1+".igenics"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                if not ioutils.checkLog4Completion("Done!", outFile+".error.log", my_env, log_file):
                    raise CraigUndefined("find_nocov_regions failed")

                shell_cmd = "getSubContigsFromIgenics.pl --cids-have-dots < "+rsqfn1+".igenics | filterRepeatedGenes.pl > "+rsqfn1+".subcontigs"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "moveLocsToSubContigs.pl -l -ctg "+rsqfn1+".subcontigs < "+chrannot_fn+" > "+rsqfn1+".locs"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "stripUTRs.pl < "+rsqfn1+".locs | filterRepeatedGenes.pl > "+rsqfn1+".nutr.locs"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "getSeqsFromGeneLocs.pl -ctg "+chrfasta_fn+" < "+rsqfn1+".subcontigs | cut -f1 > "+rsqfn1+".fa"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
# move denovo model files
                shell_cmd = "moveScoreFilterToSubContigs.pl -l -ctg "+rsqfn1+".subcontigs < "+rsqfn1+".chr.cov > "+rsqfn1+".cov"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "moveLocsToSubContigs.pl -l -ctg "+rsqfn1+".subcontigs < "+rsqfn2+".locs > "+rsqfn1+".junction.locs";
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "extractWeightedSignalFilter.pl -truncated -ctg "+rsqfn1+".fa -locs "+rsqfn1+".junction"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)

                if args.djob_num_nodes > 1:
                    ioutils.initializeDJob(args.djob_input+"/xsignals", my_env, log_file)
                    shell_cmd = "cp "+rsqfn1+".fa "+args.djob_input+"/xsignals_input"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file)
                    base_rsqfn1 = os.path.basename(rsqfn1)
                    fdtask = open(args.djob_input+"/xsignals_input/task.prop", "w")
                    fdtask.write("prefixFilesPath="+rsqfn1+"\ninputFilePath="+args.djob_input+"/xsignals_input/"+base_rsqfn1+".fa\nblockLen=20\nstranded="+tokens[5]+"\nnumPermutations="+str(args.num_permutations)+"\n")
                    fdtask.close()
                    fdcontroller = open(args.djob_input+"/xsignals_input/controller.prop", "w")
                    fdcontroller.write("masterdir="+args.djob_input+"/xsignals_input/master\ninputdir="+args.djob_input+"/xsignals_input\nnodedir="+args.djob_input+"/DJob_xsignals\nslotspernode=1\nsubtasksize=20\ntaskclass=DJob::DistribJobTasks::XSignalsTask\nnodeclass=DJob::DistribJob::"+args.djob_node_class+"\n")
                    fdcontroller.close()
                    shell_cmd = "distribjob --propFile "+args.djob_input+"/xsignals_input/controller.prop --numNodes "+str(args.djob_num_nodes)
                    subprocess.call(shell_cmd, env = my_env, stdout = log_file, stderr = log_file, shell=True)
                    if not ioutils.checkLog4Completion("Done", outFile+".error.log", my_env, log_file):
                        raise CraigUndefined("distribjob for finding xsignals failed")
                else:
                    shell_cmd = "get_covxsignals "+("--stranded" if tokens[5] == 'yes' else "")+" --prefix-evidence="+rsqfn1+" --num-permutations="+str(args.num_permutations)+" --block-len=20 "+rsqfn1+".fa"
                    ioutils.log_and_exec(shell_cmd, my_env, log_file) 
#                    print "hola"

# move utr_only model files
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
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                if not ioutils.checkLog4Completion("Done", outFile+".error.log", my_env, log_file):
                    raise CraigUndefined("joinSignalFilters.pl has found two different signals for one signle position in input")
                shell_cmd = "ln -s "+rsqfn1+".cov "+rsqfn1+".utr_only.cov"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "ln -s "+rsqfn1+".xsignal "+rsqfn1+".utr_only.xsignal"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                shell_cmd = "ln -s "+rsqfn1+".xsigscores "+rsqfn1+".utr_only.xsigscores"
                ioutils.log_and_exec(shell_cmd, my_env, log_file)
                print "hola"
            else:
                raise CraigUndefined("Unrecognized source type")
    
    if not len(prefix_files) and args.model == 'craig':
        # compute ab initio gene split of the annotation
        prefix_files = outFile
    elif not len(prefix_files):
        raise CraigUndefined("empty preconfiguration file or --model option undefined")

    base_prefix_files = os.path.basename(prefix_files)
    parmodel_prefix = os.environ['CRAIG_HOME']+"/models/"+base_prefix_files

    # obtaining a training set
    if args.model == 'ngscraig':
        fdnaive = open(parmodel_prefix+".naive.patch", "w")
        fdnaive.write("__RESOURCE_SECTION__\nResource\tev-signals-alpha\tSigma 5 -STDA -STDA - S T D A\nResource\tcoverage\tFile EXTENDED_LDSCORE PREFIX_EVIDENCE.cov true false\nResource\tRNAseq-signals\tFile XFASTA PREFIX_EVIDENCE.junction.signal ev-signals-alpha true false\nResource\tAnnot-ev-signals\tFile XFASTA PREFIX_EVIDENCE.utr_only.junction.signal ev-signals-alpha true false\nResource\tRNAseq-edges\tFile EDGEANNOT_FASTA PREFIX_EVIDENCE.junction.signal ev-signals-alpha true false\n")
        fdnaive.write("__FILTER_SECTION__\nFilter\tNorm-Coverage\tNormalizedScoreFile<double,double> coverage 1 0\nLazyFilter\tNorm-RNASeqD3-Edges\tEdgeAligner true\nFilter\tNorm-RNASeqD3-Edges-input\tFastaxFilter<EdgeAlignment,EdgeAlignment> RNAseq-edges Norm-RNASeqD3-Edges\n")
        fdnaive.write("__FEATURE_SECTION__\nFeature\tNorm-RNASeqD3-Introns-AlignType\tAlignmentType 1 Norm-RNASeqD3-Edges\nFeature\tNorm-RNASeqD3-Introns-AlignWeight\tAlignmentWeight 1 Norm-RNASeqD3-Edges\nFeature\tLinearBin4Norm-RNASeqD3-Weights\tLinearFitBin<double> Norm-RNASeqD3-Introns-AlignWeight 6 2 0 6\nPerTagFeature\tLinearBin4xNorm-RNASeqD3-IntronWeight\tFeatureXFeature<UCHAR,double> Norm-RNASeqD3-Introns-AlignType LinearBin4Norm-RNASeqD3-Weights -> CODING*, 2 INTRON_F INTRON_B, 2 INTERGENIC_F INTERGENIC_B\nFeature\tLinearBin4-Norm-Coverage\tLinearFitBin<double> null 6 2 0 6\nFeature\tNorm-Segment-Coverage\tScoringSegment<double> 3 Norm-Coverage 1 50 true\nPerTagFeature\tLinearBin4xNorm-Coverage\tFeatureOFeature<double> Norm-Segment-Coverage LinearBin4-Norm-Coverage -> CODING*, 2 INTRON_F INTRON_B, 2 INTERGENIC_F INTERGENIC_B\n")
        fdnaive.write("__PARAM_SECTION__\nOverwrite\tStart-Background\t0\t0\t5.00256\nOverwrite\tStop-Background\t0\t0\t5.0084252\nLinearBin4xNorm-RNASeqD3-IntronWeight\t4\n0\t0\n6\t0\n0\n0\n-100\n-200\n-300\n-500\n0\t0\n0\t0\nLinearBin4xNorm-RNASeqD3-IntronWeight_1	4\n1\t1\n0\t-500\n6\t0\n0\n0\n-100\n-200\n-300\n-500\n0\t0\n6\t0\n0\n10\n400\n600\n800\n1000\nLinearBin4xNorm-RNASeqD3-IntronWeight_2\t4\n0\t0\n0\t0\n6\t0\n0\n-50\n-100\n-500\n-750\n-1000\n0\t0\nLinearBin4xNorm-Coverage\t1\n6\t0\n-10\n0\n0\n0\n0\n100\nLinearBin4xNorm-Coverage_1\t1\n6\t0\n40\n20\n0\n-5\n-10\n-20\nLinearBin4xNorm-Coverage_2\t1\n6\t0\n0\n0\n0\n0\n0\n0\n")
        fdnaive.close()
        shell_cmd = "editCRAIGParamFile.pl "+parmodel_prefix+".naive.patch < "+ab_model+" > "+parmodel_prefix+".naive.params"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        shell_cmd = "craigPredict --prefix-evidence="+prefix_files+" --format=locs --DONOR-resource=RNAseq-signals --ACCEPTOR-resource=RNAseq-signals "+parmodel_prefix+".naive.params "+prefix_files+".fa > "+prefix_files+".1.locs"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        shell_cmd = "craigPredict --prefix-evidence="+prefix_files+" --format=locs --START-resource=Annot-ev-signals --STOP-resource=Annot-ev-signals --DONOR-resource=RNAseq-signals --ACCEPTOR-resource=RNAseq-signals "+parmodel_prefix+".naive.params "+prefix_files+".fa > "+prefix_files+".2.locs"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        shell_cmd = "cat "+prefix_files+".nutr.locs | grep -s -v \"<\" | grep -s -v -P \"\S>\" | filterOlapGenes.pl -v -o "+prefix_files+".2.locs | filterDiscordantIntrons.pl "+prefix_files+".2.locs >"+prefix_files+".3.locs"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

        shell_cmd = "cat "+prefix_files+".3.locs | filterDiscordantTIS.pl "+prefix_files+".1.locs > "+prefix_files+".4.locs"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

        if not args.use_non_expressed_genes:
            shell_cmd = "cat "+prefix_files+".4.locs  | filterOlapGenes.pl -s --min-weight 2 -v -o "+prefix_files+".junction.locs > "+prefix_files+".5.locs"
        else:
            shell_cmd = "mv "+prefix_files+".4.locs "+prefix_files+".5.locs"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)

        shell_cmd ="removeGenes.sh "+prefix_files+".nutr.locs "+prefix_files+".5.locs > "+prefix_files+".naivPrCorr.noAS.locs"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)
        shell_cmd = "biointers.py "+prefix_files+".fa "+prefix_files+".naivPrCorr.noAS.locs -by ic > "+prefix_files+".naivPrCorr.noAS.fa"
        ioutils.log_and_exec(shell_cmd, my_env, log_file)


        ioutils.annotateUTRs(prefix_files, "naivPrCorr.noAS", 
                             tokens[5], my_env, log_file)
        prefix_trainfiles = prefix_files+".naivPrCorr.noAS.utr"
    
    elif args.model == 'ecraig':
        prefix_trainfiles = prefix_files
    elif args.model == 'craig':
        prefix_trainfiles = prefix_files
    else:
        raise CraigUndefined("--model option undefined")

    ioutils.genCrossValRuns(prefix_trainfiles, model_name, 
                            args.out_dir, my_env, log_file)

except CraigUndefined, msg:
    print msg
