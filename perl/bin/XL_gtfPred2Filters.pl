#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Utils;
use Getopt::Long;
$| = 1;
my $cpath = "";
my $species = "";
my $usage = "usage : $0 -species species -cpath cpath < locations\n";

&GetOptions(
            "-species=s" => \$species,
            "-cpath=s" => \$cpath,
            );

(length($cpath) > 0 && length($species) > 0) || die $usage;

my $cmdOpts = "  ";
my $prg = "$ENV{CRAIG_HOME}/perl/bin/getSeqsFromLocs.pl";
my %feats = %{&Utils::groupByPattern(\*STDIN, '^>\S+\s+(\d<|<)?([^\[,\],;]+)_\d+_\d+')};

foreach my $cid (keys %feats) {
    my $ctg = $cpath."/$cid.fa";
    if(!-e $ctg) {
        $ctg = $cpath."/$cid";
    }    
    if(!-e $ctg) {
        warn "$ctg cannot be found\n";
        return;
    }
    my $myCmdOpts = $cmdOpts." -ctg $ctg ";
    Utils::execProgram($prg, $feats{$cid}, $myCmdOpts, \*STDOUT);
}



#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
$usage = "$0 species in_file in_contig out_output [in_subcontig in_augmented_contig in_augmented_marks]\n";
$species = shift || die $usage;
$file = shift || die $usage;
$contig = shift || die $usage;
$output = shift || die $usage;
$subcontig = shift;
$augmentedcontig = shift ;
$augmentedmarks = shift ;

#`estff32GTF.pl < $file.gff3 > $file.gtf`;
`gff2Locs.pl -ctg $contig < $file.gtf | grep ">" | trimStopCodon.pl -ctg $contig | grep ">" > $file.locs`;
#`genscan2Locs.pl -ctg $contig < $file.genscan | grep ">" | trimStopCodon.pl -ctg $contig | grep ">" > $file.locs`;
if(defined($subcontig)) {
    `moveLocsToSplitContigs.pl $subcontig < $file.locs > $output.pre.locs`;
}
else {
    $subcontig = $contig;
    `cp $file.locs $output.pre.locs`;
}
if(defined($augmentedcontig)) {
    `moveLocsToAugmentedContigs.pl $augmentedmarks < $output.pre.locs > $output.locs`;
}
else {
    $augmentedcontig = $subcontig;
    `cp $output.pre.locs $output.locs`;
}
`extractSignalFilter.pl -ctg $augmentedcontig < $output.locs > $output.signals`;
`filterRepeatedEsts.pl -ctg $augmentedcontig < $output.locs > $output.nr.locs`;
`$ENV{CRAIG_HOME}/bin/genes2tags --species=$species $augmentedcontig $output.nr.locs > $output.nr.ge`;
#`GEgene2est.pl $output.nr`;
`extractContentFilter.pl $augmentedcontig $output.nr.ge $output.states`;
