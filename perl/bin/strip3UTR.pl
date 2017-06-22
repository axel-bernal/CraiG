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
    undef(@{$gene->{'3UTRExons'}});
    if(scalar(@{$gene->{'Exons'}})) {
        $gene->dump(\*STDOUT, undef, 1);
    }
}

