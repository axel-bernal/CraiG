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
my $gene;
my $rgene;
my $minOlap = 1;
my $usage = "usage : $0 -min minOlapSize < locations\n";

&GetOptions(
            "min=i" => \$minOlap,
	    );

my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);

my $olapStrand = [0,0,0];
my %olapLengths = ();

if(!scalar(@genes)) {
    exit;
}


my $length;
my $lastOlap = 0;
for(my $i = 1; $i < scalar(@genes); $i++) {
    $gene = $genes[$i - 1];
    $rgene = $genes[$i];

    my $b = $rgene->compareTo($gene);

    my $olap = $minOlap - 1;
    if($gene->contigId() eq $rgene->contigId()) {
#	$olap = $gene->{'Exons'}->[-1]->[2] - $rgene->{'Exons'}->[0]->[1];
	$olap = $gene->highestCoord() - $rgene->lowestCoord() + 1;
    }

    if($minOlap <= $olap) {
	print "overlap $gene->{'Id'} $rgene->{'Id'}\n";
	if($gene->{'Strand'} ne $rgene->{'Strand'}) {
	    $olapStrand->[1]++;
	}
	else {
	    $olapStrand->[0]++;
	}
	if($lastOlap) {
	    $olapStrand->[2]++;
	}
	$olapLengths{$olap}++;            
	$lastOlap = 1;
    }
    else {
	$lastOlap = 0;
    }
}

print "same strand = ", $olapStrand->[0], "\n";
print "diff strand = ", $olapStrand->[1], "\n";
print "more than 3 genes = ", $olapStrand->[2], "\n";
print "overlap lengths\n";
#foreach my $l (sort {$a <=> $b;} keys %olapLengths) {
#    print "$l\t$olapLengths{$l}\n";
#}
