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

my %genes = % {&GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};

for my $id (keys %genes) {
    $gene = $genes{$id};
    my $exons = $gene->{'Exons'};
    $exons->[0]->[1] = $exons->[0]->[2];
#    print "->\"$gene->{'Id'} $gene->{'gId'} $gene->{'tId'}\"\n";
    if(scalar(@{$exons})) {
        $gene->dump(\*STDOUT, undef, 1);
    }
}

