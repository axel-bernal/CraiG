#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);


my $usage = "usage : $0 < locations\n";

my $gene;
my $dist;

my %genes = % {&GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};

for my $id (keys %genes) {
    $gene = $genes{$id};
#    print "->\"$gene->{'Id'} $gene->{'gId'} $gene->{'tId'}\"\n";
    my $exons = $gene->{'Exons'};
    my $fiveU = $gene->{'5UTRExons'};
    my $threeU = $gene->{'3UTRExons'};
    if(scalar(@{$fiveU}) > 0) {
	print  ">$gene->{'Id'}.5UTR\t";
	$gene->displayLoc(\*STDOUT, $gene->{'5UTRExons'});
	print "[";
	$dist = abs($fiveU->[-1]->[2] - $exons->[0]->[1]);
	if($dist > 1) {
	    print $exons->[0]->[0]."_".$exons->[0]->[1]."_".$exons->[0]->[1];
	    if(scalar(@$exons) > 1) {
		print  ">";
	    }
	}
	print "\n";
    }
    
    if(scalar(@{$threeU}) > 0) {
	print ">$gene->{'Id'}.3UTR\t";
	$dist = abs($threeU->[0]->[1] - $exons->[-1]->[2]);
	if($dist > 1) {
	    if(scalar(@$exons) > 1) {
		print  "<";
	    }
	    print $exons->[-1]->[0]."_".$exons->[-1]->[2]."_".$exons->[-1]->[2];
	}
	print "]";
	$gene->displayLoc(\*STDOUT, $gene->{'3UTRExons'});
	print "\n";
    }

    if(scalar(@{$exons})) {
	undef(@{$gene->{'5UTRExons'}});
	undef(@{$gene->{'3UTRExons'}});
	$gene->{'Id'} = $gene->{'Id'}.".CDS";
        $gene->dump(\*STDOUT, undef, 1);
    }
}
    
