#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use GenUtils;
use Getopt::Long;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 genes < contigLengths\n";

my $locs = shift || die $usage;
open(LOCS, $locs) || die "Cannot open $locs for reading\n";
my %genes = % {&GenUtils::readLocsFormat(\*LOCS)};

my %clens = ();
my @line = ();
while(<>) {
    @line = split(/\s+/);
    $clens{$line[1]} = $line[0];
}

foreach my $id (keys %genes) {
    my $gene = $genes{$id};
    my $exons = $gene->{'Exons'};
    my $strand = $gene->{'Strand'};
    if(!scalar(@{$exons})) {
	next;
    }
    my $lbegin = $exons->[0]->[1] - 100;
    my $lend = $exons->[-1]->[2] + 100;
    if($strand == STRAND_COMP) {
	$lbegin = $clens{$exons->[0]->[0]} - $exons->[0]->[1] + 1 - 100;
	$lend = $clens{$exons->[0]->[0]} - $exons->[-1]->[2] + 1 + 100;
    }
    print "$id\t$exons->[0]->[0]:$lbegin..$lend\t$strand\n";
}
