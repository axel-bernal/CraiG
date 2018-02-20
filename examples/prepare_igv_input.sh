PREFIX_OUTPUT_DIR=$1
# The prefix of the craigPreprocess output directories, i.e., without the .preproc or .model suffix
DATASET=$2
SAMPLE=$3

# call with DATASET=hehl SAMPLE=tqz
ID="$DATASET.$SAMPLE"

# Some definitions
FASTA_DIR="/Users/abernal/Data/roos8core/New"
IGVTOOLS_DIR="/Users/abernal/Applications/IGVTools"
PYTHONBIN="$CRAIG_HOME/python/bin/"
PERLBIN="$CRAIG_HOME/perl/bin/"
PREFIX_OUTPUT_DIR=$(realpath $PREFIX_OUTPUT_DIR)
PREPROC_DIR=$PREFIX_OUTPUT_DIR".preproc"
MODEL_DIR=$PREFIX_OUTPUT_DIR".model"
cd $PREPROC_DIR

#printf "Generating junction bed\n"
#$PYTHONBIN/junctions2weightedLocs.py --format bam $BAM --with-bed tgondii-rna.$ID.chr.junction.bed junction.chr.locs;

printf "Generating junction igv bed\n"
egrep -v "\t1\t-" tgondii-rna.$ID.chr.junction.bed | egrep -v "\t1\t+"  > junctions2.bed;
echo "track useScore=1" | cat - junctions2.bed > junctions.bed; rm junctions2.bed

#printf "Generating coverage\n"
#$PYTHONBIN/coverage2filter.py --format bam --len-file $FASTA_DIR/topLevelGenomicSeqs.lengths $BAM --orientation $ORIENTATION > tgondii-rna.$ID.chr.cov

printf "Generating coverage wig files\n"
$PERLBIN/cov2bed.pl $FASTA_DIR/topLevelGenomicSeqs.lengths2 depth < tgondii-rna.$ID.chr.cov
mv depth_plus.bed  depth_plus.wig;
mv depth_minus.bed  depth_minus.wig;

if [ -d $MODEL_DIR ]; then
    cd $MODEL_DIR
    $PERLBIN/locs2GTF.pl $FASTA_DIR/topLevelGenomicSeqs.fa < $MODEL_DIR/tgondi-rna.$ID.denovo.chr.locs > $MODEL_DIR/tgondi-rna.$ID.denovo.chr.gtf;
    $PERLBIN/locs2GTF.pl $FASTA_DIR/topLevelGenomicSeqs.fa < $MODEL_DIR/tgondi-rna.$ID.utr_only.chr.locs > $MODEL_DIR/tgondi-rna.$ID.utr_only.chr.gtf;
fi

printf "Building indexes for IGV browsing\n"
cd $IGVTOOLS_DIR
./igvtools   index $PREPROC_DIR/junctions.bed
./igvtools toTDF  $PREPROC_DIR/depth_plus.wig $PREPROC_DIR/depth_plus.tdf $FASTA_DIR/topLevelGenomicSeqs.fa
if [ -s $PREPROC_DIR/depth_minus.wig ]
then
    ./igvtools toTDF  $PREPROC_DIR/depth_minus.wig $PREPROC_DIR/depth_minus.tdf $FASTA_DIR/topLevelGenomicSeqs.fa
fi
