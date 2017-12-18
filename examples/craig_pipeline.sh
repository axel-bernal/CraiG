# Example of how to run craig's pipeline for Hehl's day7 RNA-Seq dataset
# This process can be made generic in a rather straightforward way.
# 1. preconfiguration file: hehl_day7.preconf which would be something like this:
# 
# evid_type(rnaseq, gpred)  type(bam|rum|gsnap|epdb*) dataset_id sample_id filename orientation(NN|N|F|R|FR|RF|FF|RR) wrapper_in
# rnaseq  bam     hehl    day7    /Users/abernal/Documents/Miscelanous/roos8core/Hehl/analyze_day7/results_sorted.bam     FR      cat
#
# all separated by tabs
#
# 2. annotation tgonME49.gtf
# 3. contigs in fasta format topLevelGenomicSeqs.fa
#

# Preprocessing the input. Need regtools in PATH
craigPreprocess.py --pre-config hehl_day7.preconf --out-dir /Users/abernal/Documents/Miscelanous/roos8core/New/hehl_day7.preproc --annot-fmt gtf --transcript-tag exon --cds-tag CDS tgondii tgonME49.gtf topLevelGenomicSeqs.fa --gc-classes=100 --model ngscraig --config config

# You can run craigReannotated directly. This will generate *.model directory
craigReannotate.py --force-train --output-model-scores --utr-only-model --prefix-output-dir hehl_day7 --model-utrs --pre-config hehl_day7.preconf  --annot-fmt gtf --transcript-tag exon --cds-tag CDS tgondii tgonME49.gtf topLevelGenomicSeqs.fa --gc-classes=100 --model ngscraig

# Run postprocessing, i.e., merging annotations made from different datasets. Run this program 
# after generating predictions from multiple datasets
craigPostprocess.py --out-dir postprocessing --list-prefixes hehl_day7.model/tgondii-rna.hehl.day7,hehl_day5.model/tgondii-rna.hehl.day5 --use-model-scores
