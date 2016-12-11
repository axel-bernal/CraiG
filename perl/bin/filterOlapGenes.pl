#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $rlocs = "";
my $reportOlap = 0;
my $sameStrand = 0;
my $oppStrand = 0;
my $minOlap = 0;
my $minWeight = 0;
my $usage = "usage : $0 --min-weight [MIN_W] [-s] [-v] -min minOlapSize -o r_locations < locations\n\t\tremoves from locations, all those genes that overlap any gene\n\t\tin r_locations\n-v\t\toutputs the overlapping set of genes instead\n-s\t\tonly reports events happening in the same strand\n-r\t\tonly reports events overlapping in opposite strands";

&GetOptions("-v" => \$reportOlap,
	    "-s" => \$sameStrand,
	    "--min-weight=i" => \$minWeight,
	    "-r" => \$oppStrand,
            "min=i" => \$minOlap,
            "-o=s" => \$rlocs,
	    );

length($rlocs) || die $usage;
open(RLOCS, "$rlocs") || die $usage;
my $strandDep = ($sameStrand ? 1 : ($oppStrand ? -1 : 0));

my %rgenes = % { &GenUtils::readLocsFormat(\*RLOCS, 0, 'OLAP', undef)};
my @rgenes = sort {&Gene::compareTo($a, $b);} (values %rgenes);

my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);

my %olapGenes = ();

if(!scalar(@genes) || !scalar(@rgenes)) {
    exit;
}

close(RLOCS);

my $i = 0;
my $length;
GenUtils::computeGeneOverlaps(\@genes, \@rgenes, \%olapGenes, $strandDep,
			      $minOlap, $minWeight);

for($i = 0; $i < scalar(@genes); $i++) {
    my $gene = $genes[$i];
    if($reportOlap && defined($olapGenes{$i})) {
	$gene->dump(\*STDOUT, undef, 1);	
    }
    elsif(!$reportOlap && !defined($olapGenes{$i})) {
	$gene->dump(\*STDOUT, undef, 1);
    }
}    

#my $reportGenes = \@genes;

#@if($reportOlap) {
 #   $reportGenes = \@olapGenes;
#}

#foreach $gene (@$reportGenes) {
#    $gene->dump(\*STDOUT, undef, 1);
#}
