#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < locations\n";

my $genes = &GenUtils::readWeightedFastaFormat(\*STDIN);
foreach my $id (keys %$genes) {
    my $exons = $genes->{$id}->{'Exons'};
    my $numExons = scalar(@{$exons});
    if(!$numExons) {
	next;
    }
    for(my $i = 0; $i < $numExons/2; $i++) {
	my $tmpExon = $exons->[$i];
	$exons->[$i] = $exons->[-($i+1)];
	$exons->[-($i+1)] = $tmpExon;
    }
	
    $genes->{$id}->displayHeader(\*STDOUT, undef);
    print "\t$exons->[0]->[3]";
    for(my $j = 1; $j < $numExons; $j++) {
	print ";", $exons->[$j]->[3];
    }
    print "\n";
}
