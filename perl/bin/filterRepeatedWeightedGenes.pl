#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $gene;
my $last_gene;
my $usage = "usage : $0 < locations\n";

my $genes = &GenUtils::readWeightedFastaFormat(\*STDIN);
my @genes = sort {&Gene::compareTo($a, $b);} (values %$genes);
if(scalar(@genes) == 0) {
    exit;
}

my $j = 0;
my $length;
my $last_length;

for(my $i = 1; $i < scalar(@genes); $i++) {
    $last_gene = $genes[$j];
    $gene = $genes[$i];
    if($gene->equalTo($last_gene)) { 
        if($last_gene->lowestCoord() == $gene->lowestCoord() 
           && $last_gene->highestCoord() == $gene->highestCoord()) {
            print STDERR "equal $last_gene->{'Id'} and $gene->{'Id'}\n";
	    for(my $e = 0; $e < scalar(@{$gene->{'Exons'}}); $e++) {
		my $exon = $gene->{'Exons'}->[$e];
		my $last_exon = $last_gene->{'Exons'}->[$e];
		$last_exon->[3] += $exon->[3];
	    }
	    
	    splice(@genes, $i, 1);
            $i--;
        }
    }
    else {
        $j = $i;
    }
}

foreach $gene (@genes) {
    $gene->displayHeader(\*STDOUT, undef);
    my $exons = $gene->{'Exons'};
    print "\t$exons->[0]->[3]";
    for(my $j = 1; $j < scalar(@$exons); $j++) {
	print ";", $exons->[$j]->[3];
    }
    print "\n";
}
