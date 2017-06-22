#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $as_events_file = "";
my $summary_locs = "";
my $models_dir = "";
my $annotation = "";
my $usage = "usage : $0 --as-events-file AS_EVENTS_FILE --summary-locs SUMMARY_LOCATION --annotation ANNOTATION --models-dir MODELS_DIRECTORY\n ";

&GetOptions(
    "--as-events-file=s" => \$as_events_file,
    "--summary-locs=s" => \$summary_locs,
    "--models-dir=s" => \$models_dir,
    "--annotation=s" => \$annotation,
    );

die $usage if $as_events_file eq "" || $annotation eq "" || $summary_locs eq "";
open(ASEVENTS, "$as_events_file") || die "Can't open file $as_events_file.\n".$usage;

my $line = "";
my %genes = ();

open(SMALL, ">small.locs") || die "Cannot open as_events.orig.locs for writing\n";
open(BIG, ">big.locs") || die "Cannot open as_events.orig.locs for writing\n";
while(defined($line = <ASEVENTS>)) {
    my @line = split(/\s/, $line);
    my $small_offset = -5;
    my $big_offset = -8;
    my $small_dbcount = $line[-2];
    my $big_dbcount = 1;
    $line[0] =~ />(\S+)/;
    my $small_eventstr = ">$1.$line[$small_offset+2]\t$line[$small_offset]"."\t".$line[$small_offset+1]."\n";
    my $big_eventstr = ">$1.$line[$big_offset+2]\t$line[$big_offset]"."\t".$line[$big_offset+1]."\n";
    print SMALL $small_eventstr;
    print BIG $big_eventstr;
    my @exons = split(/\;/, $line[$small_offset]);    
    print STDERR "missing introns for any single library for id $line[$small_offset+2].$line[0] $#exons $small_dbcount\n" if $#exons != $small_dbcount;
	
}

close(SMALL);
close(BIG);
close(ASEVENTS);

`stripUTRs.pl < $summary_locs | filterRepeatedGenes.pl | removeTIS.pl > summary.orig.nutr.locs`;
`stripUTRs.pl < $annotation | removeTIS.pl | filterRepeatedGenes.pl | biodiff.py summary.orig.nutr.locs - -by ll | biointers.py $summary_locs - -by ii | stripUTRs.pl > summary.nutr.locs`;
`clusterAnnotations.pl --mismatch-tis --match-cds summary.nutr.locs  > summary.nutr.clusters`;
`outputMatchingClusters.pl --preserve-ids < summary.nutr.clusters | grep ">" > summary.nutr.final.locs`;
`cat small.locs | filterContainedGenes.pl -v -o $summary_locs > junctions4small.locs`;
`cat big.locs | filterContainedGenes.pl -v -o $summary_locs > junctions4big.locs`;
`cat small.locs | filterContainedGenes.pl -r  -v -o $summary_locs > genes4small.locs`;
`cat big.locs | filterContainedGenes.pl -r  -v -o $summary_locs > genes4big.locs`;
