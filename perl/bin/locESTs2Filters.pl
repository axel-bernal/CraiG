#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
$usage = "$0 in_file in_contig out_output [in_subcontig in_augmented_contig in_augmented_marks]\n";
$file = shift || die $usage;
$contig = shift || die $usage;
$output = shift || die $usage;
$subcontig = shift;
$augmentedcontig = shift ;
$augmentedmarks = shift ;

`clusterESTs.pl < $file.locs | grep ">" | filterNonCannonSS.pl -ctg $contig -filter exon > $file.cluster.locs`;

$file = "$file.cluster";

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
#`filterRepeatedEsts.pl -ctg $augmentedcontig < $output.locs > $output.nr.locs`;
#`~/Projects/Devel_CRAIG/Source/genes2gelems -A $output.nr.locs -F $augmentedcontig -species celegans > $output.nr.ge`;
#`GEgene2est.pl $output.nr`;
#`extractContentFilter.pl $augmentedcontig $output.nr.est.ge $output.states`;
