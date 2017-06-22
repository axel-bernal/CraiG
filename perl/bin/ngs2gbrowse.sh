#$1 = coverage(*.chr.cov)
#$2 = junction(expecting format like junctions_all.bed)
#$3 = source -- 2nd col in the gff3
#$4 = genome
#$5 = sample.library.mapping_too
SRC=$3
LENGTHS=$4.lengths

if [ -s $3.$5.cov_plus.bw ]; 
then 
    echo "$3.$5.cov_plus.bw exists";
else
    if [ -n "$LENGTHS.gff3" ] && [ -s "$LENGTHS.gff3" ];
	echo "lengths already computed";
    then
	cat $LENGTHS | awk -v src=$SRC '{print $1 "\t" src "\tcontig\t1\t" $2 "\t.\t.\t.\tID=" $1 ";Name=" $1}' > $LENGTHS.gff3;
    fi
    cat $1 | cov2bed.pl $LENGTHS $3.$5.coverage.tmp; 
    cat $3.$5.coverage.tmp_plus.bed | awk '{printf("%s\t%d\t%d\t%.1f\n", $1, $2, $3, log($4))}' > $3.$5.cov_plus.bed;
    rm $3.$5.coverage.tmp_plus.bed;
    wigToBigWig.pl $3.$5.cov_plus.bed $LENGTHS $3.$5.cov_plus.bw;
    if [ -n $3.$5.coverage.tmp_minus.bed ];
    then
	cat $3.$5.coverage.tmp_minus.bed | awk '{printf("%s\t%d\t%d\t%.1f\n", $1, $2, $3, -log(-$4))}' > $3.$5.cov_minus.bed;
	rm $3.$5.coverage.tmp_minus.bed
	wigToBigWig.pl $3.$5.cov_minus.bed $LENGTHS $3.$5.cov_minus.bw;
    fi
fi

if [ -s $3.$5.junctions.gff3 ];
then
    echo "$3.$5.junctions.gff3 exists";
else
    FILE=$2
    if [ "${FILE##*.}" = "locs" ];
    then
	cat $2 | listIntrons.pl --include-ids | perl -ne  'if($_ =~ />(\S+)\s+(\S+)_(\d+)_(\d+)\s+(\d+)/) { print $2,"\t", '$3', "\tremark\t", $4 >= $3 ? $3 : $4, "\t", $4 >= $3 ? $4 : $3, "\t", $5, "\t", "\t+\t.\tName=", $5, "\n";}' | cat $LENGTHS.gff3 - | grep -v -w "Name=1" > $3.$5.junctions.gff3
    elif [ "${FILE##*.}" = "bed" ];
    then
	cat $2 | awk '{if($4 >= 1) {print $1 "\t" '$3' "\tremark\t" $2+51 "\t" $3-50 "\t" $4 "\t+\t.\tName=" $4}}' | sed 1,1d | cat $LENGTHS.gff3 - | grep -v -w "Name=1" > $3.$5.junctions.gff3
    fi
fi

bp_seqfeature_load.pl -a DBI::SQLite -c -f -d $3.$5.junctions.sqlite $3.$5.junctions.gff3;

if [ -d /var/lib/gbrowse2/databases/$4 ];
then
    echo "directory /var/lib/gbrowse2/databases/$4 exists!";
else
    sudo mkdir /var/lib/gbrowse2/databases/$4;
fi

sudo mv $3.$5.junctions.sqlite /var/lib/gbrowse2/databases/$4;
sudo mv $3.$5.cov_plus.bw /var/lib/gbrowse2/databases/$4/; 
if [ -e $3.$5.cov_minus.bw ];
then
    sudo mv $3.$5.cov_minus.bw /var/lib/gbrowse2/databases/$4/; 
fi
sudo mv $3.$5.junctions.sqlite /var/lib/gbrowse2/databases/$4;
sudo  chown -R apache  /var/lib/gbrowse2/databases/$4;
sudo chgrp -R apache  /var/lib/gbrowse2/databases/$4;
