#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $rlocs = "";
my $reportContainment = 0;
my $reportRLocs = 0;
my $minWeight = 0;
my $usage = "usage : $0 --min-weight [MIN_W] [-r][ -v] -o r_locations < locations\n\t\tremoves from locations, all those genes that are contained by any gene\n\t\tin r_locations\n-v\t\toutputs the contained set of genes instead\n-r\t\tOutputs the containing genes in r_locations instead\nThe containment relationship has to occur in the same strand\n";

&GetOptions("-v" => \$reportContainment,
	    "--min-weight=i" => \$minWeight,
	    "-r" => \$reportRLocs,
            "-o=s" => \$rlocs,
	    );

length($rlocs) || die $usage;

open(RLOCS, "$rlocs") || die $usage;

my %rgenes = % { &GenUtils::readLocsFormat(\*RLOCS, 0, 'OLAP', undef)};
my @rgenes = sort {&Gene::compareTo($a, $b);} (values %rgenes);

my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);

close(RLOCS);

my %contigs = ();
foreach my $id (keys %genes) {
    my $cid = $genes{$id}->{'Exons'}->[0]->[0];
    $contigs{$cid} = [[],[]] if !defined($contigs{$cid});
    push(@{$contigs{$cid}->[0]}, $genes{$id});
}

foreach my $id (keys %rgenes) {
    my $cid = $rgenes{$id}->{'Exons'}->[0]->[0];
    $contigs{$cid} = [[],[]] if !defined($contigs{$cid});
    push(@{$contigs{$cid}->[1]}, $rgenes{$id});
}

for my $cid (keys %contigs) {
    my @genes = sort {&Gene::compareTo($a, $b);} ( @{ $contigs{$cid}->[0]});
    my @rgenes = sort {&Gene::compareTo($a, $b);} ( @{ $contigs{$cid}->[1]});
    my %containment = ();
    
    next if !scalar(@genes) || !scalar(@rgenes);
    
    my $i = 0;
    my $length;
    GenUtils::computeGeneContainment(\@genes, \@rgenes, \%containment, $minWeight, $reportRLocs);
    
    my $ogenes = \@genes;
    $ogenes = \@rgenes if $reportRLocs;
    
    for($i = 0; $i <= $#$ogenes; $i++) {
	my $gene = $ogenes->[$i];
	if($reportContainment && defined($containment{$i})) {
	    $gene->dump(\*STDOUT, undef, 1);	
	}
	elsif(!$reportContainment && !defined($containment{$i})) {
	    $gene->dump(\*STDOUT, undef, 1);
	}
    }    
}

#my $reportGenes = \@genes;

#@if($reportOlap) {
 #   $reportGenes = \@olapGenes;
#}

#foreach $gene (@$reportGenes) {
#    $gene->dump(\*STDOUT, undef, 1);
#}
