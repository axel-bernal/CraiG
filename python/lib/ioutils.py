import sys, re, subprocess, logging, os

class CraigUndefined(Exception):
    pass

def clusterAnnotations(prefix, annotation, files, use_gmscore, menv, log_file):
    shell_cmd = "clusterAnnotations.pl "+(" --match-cds "+annotation+" " if annotation != "" else " ")+" ".join(files)+" > "+prefix+".clusters"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "outputMatchingClusters.pl --report-gtf "+(" --match-cds "+annotation+" " if annotation != "" else " ")+(" --use-scores " if use_gmscore else " ")+" < "+prefix+".clusters > "+prefix+".matching.clusters"
    log_and_exec(shell_cmd, menv, log_file)

    if annotation != "":
        shell_cmd = "grep -B 1 \"//\" "+prefix+".matching.clusters | grep \">\" > "+prefix+".matching.locs"
        log_and_exec(shell_cmd, menv, log_file)
        shell_cmd = "tail -1 "+prefix+".matching.clusters >> "+prefix+".matching.locs"
        log_and_exec(shell_cmd, menv, log_file)
        shell_cmd = "grep -v \"//\" "+prefix+".matching.clusters | grep -v \">\" > "+prefix+".matching.gtf"
        log_and_exec(shell_cmd, menv, log_file)

        shell_cmd = "grep -A 1 \"//\" "+prefix+".matching.clusters | grep \">\" | biodiff.py "+annotation+" - -by ll > "+prefix+".nomatching.locs"
        log_and_exec(shell_cmd, menv, log_file)
    else:
        shell_cmd = "grep \">\" "+prefix+".matching.clusters > "+prefix+".matching.locs"
        log_and_exec(shell_cmd, menv, log_file)
        shell_cmd = "grep -v \"//\" "+prefix+".matching.clusters | grep -v \">\" > "+prefix+".matching.gtf"
        log_and_exec(shell_cmd, menv, log_file)


def startCraigPredictDJob(djob_input, fa_input, model_params, prefix_data, 
                          djob_num_nodes, djob_node_class,
                          add_opts, output_file, menv, log_file):
    initializeDJob(djob_input, "craigpred", menv, log_file)
    shell_cmd = "cp "+fa_input+" "+djob_input+"/craigpred_input"
    log_and_exec(shell_cmd, menv, log_file)
    fa_input_name = os.path.basename(fa_input)
    fdtask = open(djob_input+"/craigpred_input/task.prop", "w")
    fdtask.write("paramsFile="+model_params+"\nprefixFilesPath="+prefix_data+"\ninputFilePath="+djob_input+"/craigpred_input/"+fa_input_name+"\n"+add_opts)
    fdtask.close()
    shell_cmd = "contigLength.pl -ctg "+fa_input
    num_contigs = len((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = menv, stderr = log_file).communicate()[0]).split("\n"))
    size_contigs = os.path.getsize(fa_input)
    subtask_size = int((200000*num_contigs)/size_contigs)
    subtask_size = 1 if subtask_size < 1 else subtask_size
    subtask_size = 30 if subtask_size > 30 else subtask_size       
    fdcontroller = open(djob_input+"/craigpred_input/controller.prop", "w")
    fdcontroller.write("masterdir="+djob_input+"/craigpred_input/master\ninputdir="+djob_input+"/craigpred_input\nnodedir="+djob_input+"/craigpred_DJob\nslotspernode=1\nsubtasksize="+str(subtask_size)+"\ntaskclass=DJob::DistribJobTasks::CRAIGPredTask\nnodeclass=DJob::DistribJob::"+djob_node_class+"\n")
    fdcontroller.close()
    shell_cmd = "contigLength.pl -ctg "+fa_input+" | sort -nr | head -1"
    max_contig_len = int((subprocess.Popen(shell_cmd, stdout = subprocess.PIPE, env = menv, stderr = log_file, shell = True).communicate()[0]).split()[0])
    num_gbs = 1 if max_contig_len < 6E+6 else 4
    shell_cmd = "distribjob --propFile "+djob_input+"/craigpred_input/controller.prop --numNodes "+str(djob_num_nodes)
    shell_cmd += (" --mpn "+str(num_gbs) if num_gbs > 1.8 else "")
    logging.info('Executing \"'+shell_cmd+'\"')
    subprocess.call(shell_cmd, env = menv, stdout = log_file, stderr = log_file, shell=True)
    if not checkLog4Completion("Done", menv, log_file):
        raise CraigUndefined("distribjob for CRAIGPredTask ("+prefix_data+") failed")
    
    shell_cmd = "cp "+djob_input+"/craigpred_input/master/mainresult/output0.locs "+output_file
    log_and_exec(shell_cmd, menv, log_file)


def startXSignalsDJob(djob_input, fa_input, prefix_data, 
                      djob_num_nodes, djob_node_class,
                      add_opts,
                      menv, log_file):
    initializeDJob(djob_input, "xsignals", menv, log_file)
    shell_cmd = "cp "+fa_input+" "+djob_input+"/xsignals_input"
    log_and_exec(shell_cmd, menv, log_file)
    fa_input_name = os.path.basename(fa_input)
    fdtask = open(djob_input+"/xsignals_input/task.prop", "w")
    fdtask.write("prefixFilesPath="+prefix_data+"\ninputFilePath="+djob_input+"/xsignals_input/"+fa_input_name+"\nchunkSize=30000\n"+add_opts)
    fdtask.close()
    shell_cmd = "contigLength.pl -ctg "+fa_input
    num_contigs = len((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = menv, stderr = log_file).communicate()[0]).split("\n"))
    size_contigs = os.path.getsize(fa_input)
    subtask_size = int((200000*num_contigs)/size_contigs)
    subtask_size = 1 if subtask_size < 1 else subtask_size
    subtask_size = 30 if subtask_size > 30 else subtask_size
    fdcontroller = open(djob_input+"/xsignals_input/controller.prop", "w")
    fdcontroller.write("masterdir="+djob_input+"/xsignals_input/master\ninputdir="+djob_input+"/xsignals_input\nnodedir="+djob_input+"/xsignals_DJob\nslotspernode=1\nsubtasksize="+str(subtask_size)+"\ntaskclass=DJob::DistribJobTasks::XSignalsTask\nnodeclass=DJob::DistribJob::"+djob_node_class+"\nkeepNodeForPostProcessing=yes\n")
    fdcontroller.close()
    shell_cmd = "contigLength.pl -ctg "+fa_input+" | sort -nr | head -1"
    max_contig_len = int((subprocess.Popen(shell_cmd, stdout = subprocess.PIPE, env = menv, stderr = log_file, shell = True).communicate()[0]).split()[0])
    num_gbs = 1
    if max_contig_len > 4E+7:
        num_gbs = max_contig_len/11E+6
    elif max_contig_len > 8.5E+6:
        num_gbs = max_contig_len/4E+6        
    shell_cmd = "distribjob --propFile "+djob_input+"/xsignals_input/controller.prop --numNodes "+str(djob_num_nodes)
    shell_cmd += (" --mpn "+str(num_gbs) if num_gbs > 1.8 else "")
    logging.info('Executing \"'+shell_cmd+'\"')
    subprocess.call(shell_cmd, env = menv, stdout = log_file, stderr = log_file, shell=True)
    if not checkLog4Completion("Done", menv, log_file):
        raise CraigUndefined("distribjob for XSignalTask ("+prefix_data+") failed")


def startScoreGenesDJob(djob_input, fa_input, annot_input, model_params, 
                        prefix_data, djob_num_nodes, djob_node_class,
                        add_opts, output_file,
                        menv, log_file):
    initializeDJob(djob_input, "scoregenes", menv, log_file)
    shell_cmd = "cp "+fa_input+" "+djob_input+"/scoregenes_input"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "cp "+annot_input+" "+djob_input+"/scoregenes_input"
    log_and_exec(shell_cmd, menv, log_file)
    fa_input_name = os.path.basename(fa_input)
    annot_input_name = os.path.basename(annot_input)
    fdtask = open(djob_input+"/scoregenes_input/task.prop", "w")
    fdtask.write("paramsFile="+model_params+"\nprefixEvidencePath="+prefix_data+"\nfastaFilePath="+djob_input+"/scoregenes_input/"+fa_input_name+"\nannotFilePath="+djob_input+"/scoregenes_input/"+annot_input_name+"\n"+add_opts)
    fdtask.close()
    shell_cmd = "cat "+annot_input
    num_contigs = len((subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = menv, stderr = log_file).communicate()[0]).split("\n"))
    size_contigs = os.path.getsize(fa_input)
    subtask_size = int((200000*num_contigs)/size_contigs)
    subtask_size = 1 if subtask_size < 1 else subtask_size
    subtask_size = 30 if subtask_size > 30 else subtask_size
    fdcontroller = open(djob_input+"/scoregenes_input/controller.prop", "w")
    fdcontroller.write("masterdir="+djob_input+"/scoregenes_input/master\ninputdir="+djob_input+"/scoregenes_input\nnodedir="+djob_input+"/scoregenes_DJob\nslotspernode=1\nsubtasksize="+str(subtask_size)+"\ntaskclass=DJob::DistribJobTasks::ScoreGenesTask\nnodeclass=DJob::DistribJob::"+djob_node_class+"\nkeepNodeForPostProcessing=yes\n")
    fdcontroller.close()
    shell_cmd = "contigLength.pl -ctg "+fa_input+" | sort -nr | head -1"
    max_contig_len = int((subprocess.Popen(shell_cmd, stdout = subprocess.PIPE, env = menv, stderr = log_file, shell = True).communicate()[0]).split()[0])
    num_gbs = 2 if max_contig_len < 6E+6 else 8
    shell_cmd = "distribjob --propFile "+djob_input+"/scoregenes_input/controller.prop --numNodes "+str(djob_num_nodes)
    shell_cmd += (" --mpn "+str(num_gbs) if num_gbs > 1.8 else "")
    logging.info('Executing \"'+shell_cmd+'\"')
    subprocess.call(shell_cmd, env = menv, stdout = log_file, stderr = log_file, shell=True)
    if not checkLog4Completion("Done", menv, log_file):
        raise CraigUndefined("distribjob for ScoreGenesTask ("+prefix_data+") failed")
    
    shell_cmd = "cp "+djob_input+"/scoregenes_input/master/mainresult/output.gms.locs "+output_file
    log_and_exec(shell_cmd, menv, log_file)

def computeSubContigFilters(prefix_ctg, prefix_subctg, 
                            djob_input,
                            djob_num_nodes, djob_node_class,
                            num_permutations, stranded,
                            menv, log_file) :
    shell_cmd = "moveScoreFilterToSubContigs.pl -l -ctg "+prefix_subctg+".subcontigs < "+prefix_ctg+".cov > "+prefix_subctg+".cov"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "moveLocsToSubContigs.pl -l -ctg "+prefix_subctg+".subcontigs < "+prefix_ctg+".junction.locs > "+prefix_subctg+".junction.locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "extractWeightedSignalFilter.pl -truncated -ctg "+prefix_subctg+".fa -locs "+prefix_subctg+".junction"
    log_and_exec(shell_cmd, menv, log_file)
    
    if djob_num_nodes > 1:
        add_opts = "minCoverage=2\nblockLen=20\nstranded="+("yes" if stranded else "no")+"\nnumPermutations="+str(num_permutations)+"\n"
        startXSignalsDJob(djob_input, prefix_subctg+".fa", 
                          prefix_subctg,
                          djob_num_nodes, djob_node_class,
                          add_opts, menv, log_file)
    else:
        shell_cmd = "get_covxsignals "+("--stranded" if stranded else "")+" --prefix-evidence="+prefix_subctg+" --num-permutations="+str(num_permutations)+" --chunk-size=30000 --block-len=20 "+prefix_subctg+".fa"
        log_and_exec(shell_cmd, menv, log_file) 

 
def produceSubContigAnnotation(contigs, contig_locs,
                               prefix_output, use_subcontiglocs_as_ids, preserve_ids,
                               menv, log_file) :

    shell_cmd = "cat "+contig_locs+" | extractIGenics.pl -ctg "+contigs+" > "+prefix_output+".igenics"
    log_and_exec(shell_cmd, menv, log_file)
    # using the igenics above, find subcontigs
    shell_cmd = "getSubContigsFromIgenics.pl "+("--use-subcontiglocs-as-ids " if use_subcontiglocs_as_ids else "")+"< "+prefix_output+".igenics | filterRepeatedGenes.pl > "+prefix_output+".subcontigs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "moveLocsToSubContigs.pl "+("--preserve-ids " if preserve_ids else "")+"-l -ctg "+prefix_output+".subcontigs < "+contig_locs+" > "+prefix_output+".locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "getSeqsFromGeneLocs.pl -ctg "+contigs+" < "+prefix_output+".subcontigs | cut -f1 > "+prefix_output+".fa"
    log_and_exec(shell_cmd, menv, log_file)

def produceSubsetAnnotation(contigs, igenics, contig_locs,
                            prefix_output, menv, log_file) :
    shell_cmd = "produceGeneAnnotation.pl -igenics "+igenics+" -outC "+prefix_output+".subcontigs < "+contig_locs+" > "+prefix_output+".locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "getSeqsFromGeneLocs.pl -ctg "+contigs+" < "+prefix_output+".subcontigs | cut -f1 >"+prefix_output+".fa"
    log_and_exec(shell_cmd, menv, log_file)


def annotateUTRs(prefix, prefix_files, model_name,
                 stranded, menv, log_file):
    shell_cmd = "get_transcript_utrs "+("--stranded" if stranded else "")+" -a -c --prefix-evidence="+prefix+" "+prefix_files+".fa "+prefix_files+".locs > "+prefix_files+".utr.report"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "grep -s \">\" "+prefix_files+".utr.report | filterContigBoundaryGenes.pl -ctg "+prefix_files+".fa > "+prefix_files+".utr.1.locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "biointers.py "+prefix_files+".locs "+prefix_files+".utr.1.locs -by cc > "+prefix_files+".1.locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "cat "+prefix_files+".1.locs | filterOlapGenes.pl -o "+prefix_files+".utr.1.locs | biodiff.py "+prefix_files+".1.locs - -by cc > "+prefix_files+".2.locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "biointers.py "+prefix_files+".utr.1.locs "+prefix_files+".2.locs -by cc > "+prefix_files+".utr.2.locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "produceUTRannot.pl -sort -utr "+prefix_files+".utr.2.locs < "+prefix_files+".2.locs > "+prefix_files+".utr.3.locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "biointers.py "+prefix+".fa "+prefix_files+".utr.3.locs -by ic > "+prefix_files+".utr.3.fa"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "genes2tags --DONOR-resource=RNAseq-signals --ACCEPTOR-resource=RNAseq-signals --have-utrs --species="+model_name+" --prefix-evidence="+prefix+" "+prefix_files+".utr.3.fa "+prefix_files+".utr.3.locs 2> "+prefix_files+".utr.err > "+prefix_files+".utr.3.tag"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "grep -s -v \"START_\" < "+prefix_files+".utr.err | grep -s -v \"STOP_\" | grep -s \"Added\" | perl -ne '$_ =~ /.+\s(\S+)$/; print \"$1\n\";' | sort -u | awk '{print \">\" $1}' > "+prefix_files+".utr.3.ids"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "biodiff.py "+prefix_files+".utr.3.locs "+prefix_files+".utr.3.ids -by ci > "+prefix_files+".utr.locs"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "biodiff.py "+prefix_files+".utr.3.fa "+prefix_files+".utr.3.ids -by ii > "+prefix_files+".utr.fa"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "biodiff.py "+prefix_files+".utr.3.tag "+prefix_files+".utr.3.ids -by ii > "+prefix_files+".utr.tag"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "rm "+prefix_files+".utr.?.locs "+prefix_files+".utr.?.fa "+prefix_files+".utr.?.tag "+prefix_files+".utr.?.ids"
    log_and_exec(shell_cmd, menv, log_file)


def createSetup4Training(prefix_train, prefix_genome, 
                         prefix_xvalfiles, 
                         model_name, gc_classes,
                         out_dir, conf_file,
                         menv, log_file) :
    #creating cross validation runs
    base_xvalfiles = os.path.basename(prefix_xvalfiles)
    shell_cmd = "cd "+out_dir+" && crossValidRuns.pl -notest 12 "+base_xvalfiles+".locs "+base_xvalfiles+".fa "+base_xvalfiles+".tag"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "cd "+out_dir+" && processAnnotation.pl -s "+base_xvalfiles+".run.0/"+base_xvalfiles+".valid "+base_xvalfiles+".run.0/"+base_xvalfiles+".valid "+gc_classes
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "grep -v Broken "+log_file.name+" > "+log_file.name+"2; mv "+log_file.name+"2 "+log_file.name
    log_and_exec(shell_cmd, menv, log_file)
    # creating config file
    shell_cmd = "cd "+out_dir+" && pwd"
    absout_dir = subprocess.Popen(shell_cmd, stdout = subprocess.PIPE, env = menv, stderr = log_file, shell=True).communicate()[0].rstrip('\n')
    configh = open(absout_dir+"/"+conf_file, "w")
    configh.write("ValidationSequences "+absout_dir+"/"+base_xvalfiles+".run.0/"+base_xvalfiles+".valid.fa\n")
    configh.write("ValidationTags "+absout_dir+"/"+base_xvalfiles+".run.0/"+base_xvalfiles+".valid.tag\n")
    configh.write("Sequences "+absout_dir+"/"+base_xvalfiles+".run.0/"+base_xvalfiles+".train.fa\n")
    configh.write("Tags "+absout_dir+"/"+base_xvalfiles+".run.0/"+base_xvalfiles+".train.tag\n")
    configh.write("Name "+model_name+"\n")
    configh.write("Path "+absout_dir+"/"+base_xvalfiles+".run.0\n")
    configh.write("PrefixFiles "+absout_dir+"/"+os.path.basename(prefix_train)+"\n")
    configh.write("PrefixGenomeFiles "+absout_dir+"/"+os.path.basename(prefix_genome)+"\n")
    configh.write("PrefixXValFiles "+absout_dir+"/"+os.path.basename(prefix_xvalfiles)+"\n")
    configh.close()


def checkLog4Completion(phrase, menv, log_file) :
    shell_cmd = "tail -1 "+log_file.name
    prg_output = subprocess.Popen(shell_cmd.split(), stdout = subprocess.PIPE, env = menv, stderr = log_file).communicate()[0]
    match = re.match(r"^"+phrase, prg_output.split("\n")[0])
    return match


def log_and_exec(shell_cmd, menv, log_file) :
    logging.info('Executing \"'+shell_cmd+'\"')
    return_code = subprocess.call(shell_cmd, env = menv, stderr = log_file, shell=True)
    logging.info('Done \"'+shell_cmd+'\"!')
    return return_code

def processFasta(script_in, fmt, fn_in, fn_out, menv, log_file):
    logging.info('Processing sequence file in '+fmt+' format')
    shell_cmd = script_in+" < "+fn_in+" > "+fn_out+".0.fa"
    log_and_exec(shell_cmd, menv, log_file)

    if fmt == 'genbank':        
        shell_cmd = "genbank2fa.pl < "+fn_out+".0.fa > "+fn_out+".fa"
        log_and_exec(shell_cmd, menv, log_file)
    elif fmt == 'fasta':
        # make contigs to be one long line. Useful for parsing
        shell_cmd = "perl -ne \'chomp(); if($_ =~ /^>(\S+)\s*/) { print \"\n>$1\n\";} else { print $_; }\' < "+fn_out+".0.fa | sed 1,1d | sed -e '$a\\' > "+fn_out+".fa"
        log_and_exec(shell_cmd, menv, log_file)
    else:
        raise CraigUndefined(fmt+" is not defined")

    shell_cmd = "rm "+fn_out+".?.fa"
    log_and_exec(shell_cmd, menv, log_file)

def processAnnotation(script_in, fmt, tr_tag, cds_tag, non_uniq_ids, fn_in, fn_out, menv, log_file):
    logging.info('Processing annotation file in '+fmt+' format')
    shell_cmd = script_in+" < "+fn_in+ " > "+fn_out+".0.locs"
    log_and_exec(shell_cmd, menv, log_file)

    gtfBin = "gtf2Locs.pl "+("" if non_uniq_ids else "--uniq-id")
    
    if fmt == 'gff3':
        shell_cmd = "grep -s -P -v \"\s"+tr_tag+"\s\" "+fn_out+".0.locs | gff32GTF.pl -all | "+gtfBin+" | sort > "+fn_out+".1.locs"
        log_and_exec(shell_cmd, menv, log_file)
        if tr_tag != 'null':
            shell_cmd = "grep -s -P -v \"\s"+cds_tag+"\s\" "+fn_out+".0.locs | gff32GTF.pl -all | "+gtfBin+" > "+fn_out+".3.locs"
            log_and_exec(shell_cmd, menv, log_file)
    elif fmt == 'gtf':
        shell_cmd = "grep -s -P -v \"\s"+tr_tag+"\s\" "+fn_out+".0.locs | "+gtfBin+" | sort > "+fn_out+".1.locs"
        log_and_exec(shell_cmd, menv, log_file)
        if tr_tag != 'null':
            shell_cmd = "grep -s -P -v \"\s"+cds_tag+"\s\" "+fn_out+".0.locs | "+gtfBin+" > "+fn_out+".3.locs"
            log_and_exec(shell_cmd, menv, log_file)
    elif fmt == genbank:
        subprocess.call("genbank2Locs.pl --cds-tag "+cds_tag+" --transcript-tag "+tr_tag+" < "+fn_out+".0.locs > "+fn_out+".1.locs", env = menv, stderr = log_file, shell=True)
        log_and_exec(shell_cmd, menv, log_file)
    elif fmt == locs:
        shell_cmd = "cp "+fn_in+" "+fn_out+".1.locs"
        log_and_exec(shell_cmd, menv, log_file)
    else:
        raise CraigUndefined(fmt+" is not defined")

    if tr_tag != 'null' and (fmt == 'gff3' or fmt == 'gtf'):
        shell_cmd = "biointers.py "+fn_out+".3.locs "+fn_out+".1.locs -by ii | sort > "+fn_out+".4.locs"
        log_and_exec(shell_cmd, menv, log_file)
        shell_cmd = "produceUTRannot.pl -utr "+fn_out+".4.locs < "+fn_out+".1.locs > "+fn_out+".2.locs; mv "+fn_out+".2.locs "+fn_out+".1.locs"
        log_and_exec(shell_cmd, menv, log_file)

    shell_cmd = "filterRepeatedGenes.pl < "+fn_out+".1.locs > "+fn_out+".locs; rm "+fn_out+".?.locs"
    log_and_exec(shell_cmd, menv, log_file)
    
    fixTranslatedProduct(fn_out+".fa", fn_out+".locs", menv, log_file)


def computePreprocAuxFiles(chrfasta_fn, chrannot_fn, chrlen_fn, chrlen2_fn, chrss_fn, menv, log_file):
    shell_cmd = "contigLength.pl -ctg "+chrfasta_fn+" > "+chrlen_fn
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "perl -ne '$_ =~ /^(\S+)\s+(\S+)\s*$/; print \">$2\t$1\n\/\/\n\";' < "+chrlen_fn+" > "+chrlen2_fn
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "reportSSpairs.pl -ctg "+chrfasta_fn+" < "+chrannot_fn+" > "+chrss_fn
    log_and_exec(shell_cmd, menv, log_file)
    

def fixTranslatedProduct(chrfasta_fn, chrannot_fn, menv, log_file):
    shell_cmd = "cat "+chrannot_fn+" | fixTranslatedProduct.pl --preserve-phases -ctg "+chrfasta_fn+" | grep -s \">\" > "+chrannot_fn+".1 ; mv "+chrannot_fn+".1 "+chrannot_fn
    log_and_exec(shell_cmd, menv, log_file)

# strict_remove is true if we want to remove contigs in which at least one annotation was removed. 
# A value of false will only remove those contigs which have no annotations in chrannot_fn
def synchronizeAnnotation(strict_remove, chrfasta_fn, oldcharannot_fn, chrannot_fn, menv, log_file):
    shell_cmd = "biointers.py "+chrfasta_fn+" "+chrannot_fn+" -by ic > "+chrfasta_fn+".1"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "biointers.py "+chrannot_fn+" "+chrfasta_fn+".1 "+chrannot_fn+" -by ic > "+chrannot_fn+".1"
    log_and_exec(shell_cmd, menv, log_file)

def initializeDJob(djob_inpdir, djob_prefix, menv, log_file):
    if not os.path.exists(djob_inpdir):
        shell_cmd = "mkdir "+djob_inpdir
        log_and_exec(shell_cmd, menv, log_file)

    djob_name = djob_inpdir+"/"+djob_prefix
    shell_cmd = "rm -rf "+djob_name+"_DJob/*" if os.path.exists(djob_name+"_DJob") else "mkdir "+djob_name+"_DJob"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "rm -rf "+djob_name+"_input/*" if os.path.exists(djob_name+"_input") else "mkdir "+djob_name+"_input" 
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "mkdir "+djob_name+"_input/master"
    log_and_exec(shell_cmd, menv, log_file)

def makeNaiveRSqParamModel(naive_model, abi_model, menv, log_file):
    fdnaive = open(naive_model+".patch", "w")
    fdnaive.write("__RESOURCE_SECTION__\nResource\tev-signals-alpha\tSigma 5 -STDA -STDA - S T D A\nResource\tcoverage\tFile EXTENDED_LDSCORE PREFIX_EVIDENCE.cov true false\nResource\tRNAseq-signals\tFile XFASTA PREFIX_EVIDENCE.junction.signal ev-signals-alpha true false\nResource\tAnnot-ev-signals\tFile XFASTA PREFIX_EVIDENCE.utr_only.junction.signal ev-signals-alpha true false\nResource\tRNAseq-edges\tFile EDGEANNOT_FASTA PREFIX_EVIDENCE.junction.signal ev-signals-alpha true false\n")
    fdnaive.write("__FILTER_SECTION__\nFilter\tNorm-Coverage\tNormalizedScoreFile<double,double> coverage 1 0\nLazyFilter\tNorm-RNASeqD3-Edges\tEdgeAligner true\nFilter\tNorm-RNASeqD3-Edges-input\tFastaxFilter<EdgeAlignment,EdgeAlignment> RNAseq-edges Norm-RNASeqD3-Edges\n")
    fdnaive.write("__FEATURE_SECTION__\nFeature\tNorm-RNASeqD3-Introns-AlignType\tAlignmentType 1 Norm-RNASeqD3-Edges\nFeature\tNorm-RNASeqD3-Introns-AlignWeight\tAlignmentWeight 1 Norm-RNASeqD3-Edges\nFeature\tLinearBin4Norm-RNASeqD3-Weights\tLinearFitBin<double> Norm-RNASeqD3-Introns-AlignWeight 6 2 0 6\nPerTagFeature\tLinearBin4xNorm-RNASeqD3-IntronWeight\tFeatureXFeature<UCHAR,double> Norm-RNASeqD3-Introns-AlignType LinearBin4Norm-RNASeqD3-Weights -> CODING*, CDS_INTRONS*, 2 INTERGENIC_F INTERGENIC_B\nFeature\tLinearBin4-Norm-Coverage\tLinearFitBin<double> null 6 2 0 6\nFeature\tNorm-Segment-Coverage\tScoringSegment<double> 3 Norm-Coverage 1 50 true\nPerTagFeature\tLinearBin4xNorm-Coverage\tFeatureOFeature<double> Norm-Segment-Coverage LinearBin4-Norm-Coverage -> CODING*, CDS_INTRONS*, 2 INTERGENIC_F INTERGENIC_B\n")
    fdnaive.write("__PARAM_SECTION__\nOverwrite\tStart-Background\t0\t0\t5.00256\nOverwrite\tStop-Background\t0\t0\t5.0084252\nLinearBin4xNorm-RNASeqD3-IntronWeight\t4\n0\t0\n6\t0\n0\n0\n-100\n-200\n-300\n-500\n0\t0\n0\t0\nLinearBin4xNorm-RNASeqD3-IntronWeight_1	4\n1\t1\n0\t-500\n6\t0\n0\n0\n-100\n-200\n-300\n-500\n0\t0\n6\t0\n0\n10\n400\n600\n800\n1000\nLinearBin4xNorm-RNASeqD3-IntronWeight_2\t4\n0\t0\n0\t0\n6\t0\n0\n-50\n-100\n-500\n-750\n-1000\n0\t0\nLinearBin4xNorm-Coverage\t1\n6\t0\n-10\n0\n0\n0\n0\n100\nLinearBin4xNorm-Coverage_1\t1\n6\t0\n40\n20\n0\n-5\n-10\n-20\nLinearBin4xNorm-Coverage_2\t1\n6\t0\n0\n0\n0\n0\n0\n0\n")
    fdnaive.close()
    shell_cmd = "editCRAIGParamFile.pl "+naive_model+".patch < "+abi_model+" > "+naive_model+".tmp"
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "mv "+naive_model+".tmp "+naive_model
    log_and_exec(shell_cmd, menv, log_file)
    shell_cmd = "rm "+naive_model+".patch"
    log_and_exec(shell_cmd, menv, log_file)
        
