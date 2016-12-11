lib_LIBRARIES = libcraig.a
libcraig_a_SOURCES = Overlap.cpp Overlap.h GeneTagPrinter.h  RepeatedElem.h GeneTagPrinter.cpp Sim.cpp Exon.cpp RSqUtils.h RSqUtils.cpp GeneUtils.h GeneUtils.cpp Sim.h Exon.h Gene.cpp GeneEvaluator.cpp Transcript.cpp Transcript.h TypeDefs.cpp Gene.h Organism.cpp Organism.h TypeDefs.h GeneEvaluator.h ArgParseUtils.h

bin_PROGRAMS = craigPredict scraigPredict craigTrain scraigTrain rcraigTrain genes2tags train_cimm compute_feat_vector emit_sequence gene_extender list_unreachable_genes extract_orfs score_genes score_tags prepare_rnaseq_junctions ngscraigPredict report_rnaseq_asevents find_nocov_regions eval_prediction compute_param_stats get_covxsignals get_transcript_utrs get_transcript_utrs_report analyze_annot_vs_rnaseq aggregate_scorefilter_vals get_coverage_genes aggregate_signalfilter_vals remove_suspect_transcripts 

craigPredict_SOURCES = craigPredict.cpp
scraigPredict_SOURCES = scraigPredict.cpp
craigTrain_SOURCES = craigTrain.cpp
scraigTrain_SOURCES = scraigTrain.cpp
rcraigTrain_SOURCES = rcraigTrain.cpp
genes2tags_SOURCES = genes2tags.cpp
train_cimm_SOURCES = train_cimm.cpp
compute_feat_vector_SOURCES = compute_feat_vector.cpp
emit_sequence_SOURCES = emit_sequence.cpp
gene_extender_SOURCES = gene_extender.cpp
list_unreachable_genes_SOURCES = list_unreachable_genes.cpp
extract_orfs_SOURCES = extract_orfs.cpp
score_tags_SOURCES = score_tags.cpp
score_genes_SOURCES = score_genes.cpp
prepare_rnaseq_junctions_SOURCES = prepare_rnaseq_junctions.cpp
ngscraigPredict_SOURCES = ngscraigPredict.cpp
report_rnaseq_asevents_SOURCES = report_rnaseq_asevents.cpp
find_nocov_regions_SOURCES = find_nocov_regions.cpp
eval_prediction_SOURCES = eval_prediction.cpp
compute_param_stats_SOURCES = compute_param_stats.cpp
get_covxsignals_SOURCES = get_covxsignals.cpp
get_transcript_utrs_SOURCES = get_transcript_utrs.cpp
get_transcript_utrs_report_SOURCES = get_transcript_utrs_report.cpp
analyze_annot_vs_rnaseq_SOURCES = analyze_annot_vs_rnaseq.cpp
aggregate_scorefilter_vals_SOURCES = aggregate_scorefilter_vals.cpp
aggregate_signalfilter_vals_SOURCES = aggregate_signalfilter_vals.cpp
get_coverage_genes_SOURCES = get_coverage_genes.cpp
remove_suspect_transcripts_SOURCES = remove_suspect_transcripts.cpp
AM_LDFLAGS = @top_srcdir@/lless/liblless.la @top_srcdir@/regex/libboost_regex.a

AM_CPPFLAGS = -I@top_srcdir@/src/ -I@top_srcdir@/lless/ -I@top_srcdir@/regex -I@top_srcdir@/sparsehash/include

include_HEADERS =  Overlap.h GeneTagPrinter.h  RepeatedElem.h RSqUtils.h GeneUtils.h Sim.h Exon.h  Transcript.h Gene.h Organism.h TypeDefs.h GeneEvaluator.h ArgParseUtils.h


LDADD = libcraig.a
