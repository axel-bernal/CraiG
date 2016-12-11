#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

my $usage = "$0 locs contigs locs2keep [-soft]\nOption -soft indicates that alternative spliced transcripts that are removed do not cause the gene to be removed. locs must be a tab-separated file\n";
my $locs = shift || die $usage;
my $contigs = shift || die $usage;
my $toKeep = shift || die $usage;
my $soft = shift;

`sort $locs > $locs.$$`;
`grep ">" $contigs | cut -f2 -d ">" | cut -f1 |sort -u  >  $contigs.$$`;
if(defined($soft) && $soft eq "-soft") {
    `cut -f1 -d "|" $toKeep | sort -u | join -t "|" -v 2 - $locs.$$ | cut -f2 | cut -f2 -d "<" | cut -f1 -d "_" | sort -u > $contigs.toRemove.$$`;
}    
else {
#    `sort $toKeep | join -v 2 - $locs.$$ | cut -f2 -d " " | cut -f2 -d "<" | cut -f1 -d "_" | sort -u > $contigs.toRemove.$$`;
    `cut -f1 $toKeep | cut -f1 -d " " | cut -f2 -d ">" | joinBioLists.pl fasta $locs.$$ -v | listContigs.pl | sort -u > $contigs.toRemove.$$`;
}

`joinBioLists.pl locs $locs < $contigs.toRemove.$$ | cut -f1 | sort -u > $locs.toRemove.$$`;
#`cut -f1 $toKeep | sort | join -v2 - $locs.$$ >> $locs.toRemove.$$`
`cut -f1 $toKeep | cut -f2 -d ">"  | joinBioLists.pl fasta $locs.$$ -v >> $locs.toRemove.$$`;
`sort -u $locs.toRemove.$$ > $locs.toRemove2.$$`;
#`cp  $locs.toRemove2.$$ $locs.toRemove.$$`;
#`join -v 1 $locs.$$ $locs.toRemove.$$ > $locs.new`;
`cut -f1 $locs.toRemove2.$$ | cut -f2 -d ">" > $locs.toRemove.$$`;
`joinBioLists.pl fasta $locs.$$ -v < $locs.toRemove.$$ > $locs.new`;

`join -v 2 $contigs.toRemove.$$ $contigs.$$ | sort -u > $contigs.toKeep.$$`;
`joinBioLists.pl fasta $contigs < $contigs.toKeep.$$ > $contigs.new`;
`rm $locs.$$ $contigs.$$ $locs.toRemove.$$ $locs.toRemove2.$$ $contigs.toRemove.$$ $contigs.toKeep.$$`;
