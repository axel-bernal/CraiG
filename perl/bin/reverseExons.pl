#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < locations\n";

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
foreach my $id (keys %genes) {
    my $exons = $genes{$id}->{'Exons'};
    my $numExons = scalar(@{$exons});

    for(my $i = 0; $i < $numExons/2; $i++) {
	reverseExon($exons->[$i]);
	reverseExon($exons->[-($i+1)]);
	my $tmpExon = $exons->[$i]; 
	$exons->[$i] = $exons->[-($i+1)];
	$exons->[-($i+1)] = $tmpExon;	
    }

    $genes{$id}->dump(\*STDOUT);
}

sub reverseExon {
    my $exon = shift;
    my $tmp = $exon->[1];
    $exon->[1] = $exon->[2];
    $exon->[2] = $tmp;
}
