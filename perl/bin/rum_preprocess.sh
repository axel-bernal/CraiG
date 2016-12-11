#$1(genome in fa) = ME49.fa
#$2(annotation in gff) = TgondiiME49_ToxoDB-7.0.gff
# $3(NAME) = me49
# $4(source) = db8.0
# $5(cds keyword) = CDS|exon
BASE_DIR=$HOME/Applications/Contrib/rum/index_maker;
formatFastaFile.pl < $1 > $3.chr.fa;
ln -s $3.chr.fa $3_genome.txt;
cat $2 | grep -P "\t$5\t" | gff32bed.pl > $3_$4.txt;
perl -ne 'chomp; if($_ !~ /\#\#/) { @line = split(/\t/, $_); @starts = split(/,/, $line[3]); @ends = split(/,/, $line[4]); print "$line[1]\t$line[2]\t$starts[0]\t$ends[-1]\t", scalar(@starts), "\t$line[3]\t$line[4]\t$line[0]\n";}' < $3_$4.txt > $3_$4_cluster.txt
CURR_DIR=`pwd`
echo $CURR_DIR/"$3_$4.txt" > gene_info_files;
cp gene_info_files ~/Applications/Contrib/rum/index_maker;

cd $BASE_DIR && perl ./create_indexes_from_ucsc.pl $CURR_DIR/$3_genome.txt $CURR_DIR/$3_$4;
cd $CURR_DIR;
mv $3_$4_gene_info.txt ~/Applications/Contrib/rum/indexes
mv $3_genome_one-line-seqs.fa ~/Applications/Contrib/rum/indexes
mv $3_genes.4.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genes.3.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genes.2.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genes.1.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genes.rev.2.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genes.rev.1.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genome.4.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genome.3.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genome.2.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genome.1.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genome.rev.2.ebwt ~/Applications/Contrib/rum/indexes
mv $3_genome.rev.1.ebwt ~/Applications/Contrib/rum/indexes

#cd ~/Applications/Contrib/rum/ && nohup perl RUM_runner.pl lib/rum.config_$3 /home/abernal/Data/toxodb-7.1/$1/RNA-Seq/oocyst_mrna.fastq data/$1.oocyst 2 oocyst 2> data/$1.oocyst/rum.running.log
