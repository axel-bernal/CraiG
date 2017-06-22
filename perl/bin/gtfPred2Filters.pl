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

`gtf2Locs.pl -ctg $contig < $file.gtf 2> $file.log | grep ">" | trimStopCodon.pl $contig | grep ">" > $file.locs`;
#`genscan2Locs.pl $contig < $file.genscan | grep ">" | trimStopCodon.pl $contig | grep ">" > $file.locs`;
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
`filterRepeatedGenes.pl -ctg $augmentedcontig < $output.locs > $output.nr.locs`;
`$ENV{CRAIG_HOME}/bin/genes2tags --species=$species $augmentedcontig $output.nr.locs > $output.nr.ge`;

$output = "$output.nr";
`extractContentFilter.pl $augmentedcontig $output.ge $output.states`;
