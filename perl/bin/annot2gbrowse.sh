#$1 = annotation
#$2 = fasta
#$3 = source -- this goes in the gtf as second column
#$4 = genome -- name and version of genome
#$5 = tag for annotation -- use annot for official annotation or sample.library.mapping_tool.prediction_tool for RNASeq derived or source/author.prediction_tool for gene prediction
#$6 database for gbrowse stub
SRC=$3
LENGTHS=$4.lengths

if [ -n "$LENGTHS.gff3" ] && [ -s "$LENGTHS.gff3" ];
then
    echo "lengths already computed";
else
    contigLength.pl -ctg $2 | awk '{print $2 "\t" $1}' > $LENGTHS;
    cat $LENGTHS | awk -v src=$SRC '{print $1 "\t" src "\tcontig\t1\t" $2 "\t.\t.\t.\tID=" $1 ";Name=" $1}' > $LENGTHS.gff3;
fi

FILE=$1
if [ "${FILE##*.}" = "gff3" ];
then
    cat $FILE | sed 's/CRAIG/'$3'/g' | sed 's/5UTR/UTR/g' | sed 's/3UTR/UTR/g' > $3.$5.tmp.gff3;
elif [ "${FILE##*.}" = "gff" ]; 
then
    gff32GTF.pl -all < $FILE | gtf2Locs.pl -ctg $2 -uniq-id | fixTranslatedProduct.pl -ctg $2 | grep ">" | locs2GTF.pl $2 | gtf2gff3.pl | sed 's/CRAIG/'$3'/g' | sed 's/5UTR/UTR/g' | sed 's/3UTR/UTR/g' > $3.$5.tmp.gff3;
else   
    cat $FILE | locs2GTF.pl $2 | gtf2gff3.pl | sed 's/CRAIG/'$3'/g' | sed 's/5UTR/UTR/g' | sed 's/3UTR/UTR/g' > $3.$5.tmp.gff3;
fi

sed 1,1d $3.$5.tmp.gff3 | cat $LENGTHS.gff3 - > $3.$5.gff3;
rm $3.$5.tmp.gff3;
bp_seqfeature_load.pl -a DBI::SQLite -c -f -d $3.$5.sqlite $3.$5.gff3;

if [ -d /var/lib/gbrowse2/databases/$4 ];
then
    echo "directory /var/lib/gbrowse2/databases/$4 exists!";
else
    sudo mkdir /var/lib/gbrowse2/databases/$4;
fi

sudo mv $3.$5.sqlite /var/lib/gbrowse2/databases/$4;
sudo chown -R apache  /var/lib/gbrowse2/databases/$4;
sudo chgrp -R apache  /var/lib/gbrowse2/databases/$4;

