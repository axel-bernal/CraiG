#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);


my $usage = "usage : $0 contig < locations\n";

my $gene;
my $reportSeq = 0;
&GetOptions("-reportSeq" => \$reportSeq,);

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};

foreach my $id (keys %genes) {
    $gene = $genes{$id};
    print $gene->{'Id'}, "\t";
    print $gene->contigId(),"_", $gene->lowestCoord(), "_";
    print $gene->highestCoord(), "\t";
    print scalar(@{$gene->{'5UTRExons'}}), "\t";
    
    if(scalar(@{$gene->{'5UTRExons'}}) > 0) {
	print $gene->getFivePWeight(), "\t";
	my $length = 0;
	for(my $j = 0; $j < scalar(@{$gene->{'5UTRExons'}}); $j++) {	
	    $length += abs($gene->{'5UTRExons'}->[$j]->[2] - $gene->{'5UTRExons'}->[$j]->[1]) + 1; 
	}
	print $length,"\t";
    }	
    else {
	print "0\t0\t";
    }

    print scalar(@{$gene->{'3UTRExons'}}), "\t";
    if(scalar(@{$gene->{'3UTRExons'}}) > 0) {
	print $gene->getThreePWeight(), "\t";
	my $length = 0;
	for(my $j = 0; $j < scalar(@{$gene->{'3UTRExons'}}); $j++) {	
	    $length += abs($gene->{'3UTRExons'}->[$j]->[2] - $gene->{'3UTRExons'}->[$j]->[1]) + 1; 
	}
	print $length,"\n";
    }
    else {
	print "0\t0\n";
    }
}

