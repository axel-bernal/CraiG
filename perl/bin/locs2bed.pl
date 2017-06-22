#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use GenUtils;
use Getopt::Long;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 < genes\n";

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};

print "##name\tchrom\tstrand\texonStarts\texonEnds\n";
foreach my $id (keys %genes) {
    my $gene = $genes{$id};
    my $exons = $gene->{'Exons'};
    my $strand = ($gene->{'Strand'} == STRAND_FWD ? "+" : "-");
    if(!scalar(@{$exons})) {
	next;
    }
    print "$id\t$exons->[0]->[0]\t$strand\t";
    my $i;
    if($strand eq '+') {
	for($i = 0; $i < scalar(@$exons); $i++) {
	    print $exons->[$i]->[1] - 1;
	    if($i < scalar(@$exons) - 1) {print ",";};
	}
    }
    else {
	for($i = scalar(@$exons) - 1; $i >= 0; $i--) {
	    print $exons->[$i]->[2] - 1;
	    if($i > 0) {print ",";};
	}
    }
    print "\t";
    if($strand eq '+') {
	for($i = 0; $i < scalar(@$exons); $i++) {
	    print $exons->[$i]->[2];
	    if($i < scalar(@$exons) - 1) {print ",";};
	}
    }
    else {
	for($i = scalar(@$exons) - 1; $i >= 0; $i--) {
	    print $exons->[$i]->[1];
	    if($i > 0) {print ",";};
	}	
    }
    print "\n";
}
