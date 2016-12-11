#$1(genome in fa) = ME49.fa
#$2(annotation in gff) = TgondiiME49_ToxoDB-7.0.gff
# $3(NAME) = me49
# $4(cds keyword) = CDS|exon
formatFastaFile.pl < $1 |  fasta2files.pl -dir $3.fa;
#buffer_size=`eval ls -ltr $3.fa/* | wc | awk '{print $3}'`
ls $3.fa/* > $3.fa.list
gmap_build -d $3 -k 15 $3.fa.list;
cat $2 | grep -P "\t$4\t" | gff32GTF.pl -all | sed 's/\tCDS\t/\texon\t/g' | gtf_splicesites > $3.splicesites;
cat $3.splicesites | iit_store -o $3.splicesites.iit;
cp $3.splicesites.iit ~/Applications/Contrib/gmap/share/$3/$3.maps;
#gsnap -d  me49 -k 15 -s me49.splicesites.iit -N 1 --split-output=oocyst_gsnap -t 4 ../oocyst_mrna.fastq 2> oocyst_gsnap.output.log
