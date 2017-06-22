#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";

my $preserve_phases = 0; # recompute phases..  should be given as an option.. TODO
my $usage = "$0 [-evid] species contigs annot\n";
my $evid = 0;
my $species = shift || die $usage;

if($species eq "-evid") {
    $evid = 1;
    $species = shift || die $usage;
}

my $contigs = shift || die $usage;
my $annot = shift || die $usage;
my $ext = "gff";

(-e "$annot") || die "Can't find $annot.\n$usage";

if($annot =~ /(\S+)\.([^\.]+)/) {
    $annot = $1;
    $ext = $2;
}
else {
    die "$annot needs an extension.\n$usage";
}


if($ext eq "gtf" || $ext eq "gff") {
    `gtf2Locs.pl -ctg $contigs < $annot.$ext 2> $annot.log | grep ">" > $annot.nr.locs`;
}
else {
    system("checkAnnotSanity.pl ".($preserve_phases ? "--preserve-phases " : "")."-ctg $contigs < $annot.locs 2> $annot.log | grep ">" > $annot.nr.locs");
}

`filterRepeatedGenes.pl < $annot.nr.locs 2>> $annot.log > $annot.locs; rm $annot.nr.locs`;

if($evid) {
#    if(!-e "$annot.tag" || -z "$annot.tag") {
        `genes2tags --species=$species $contigs $annot.locs > $annot.tag 2>> $annot.log`;
#    }
#    if(!-e "$annot.state" || -z "$annot.state") {
        `extractContentFilter.pl -ctg $contigs < $annot.tag > $annot.state 2>> $annot.log`;
#    }
#    if(!-e "$annot.signal" || -z "$annot.signal") {
        `extractSignalFilter.pl -ctg $contigs < $annot.locs > $annot.signal 2>> $annot.log`;
#    }
}
