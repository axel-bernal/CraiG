#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Gene;
use GenUtils;
use Getopt::Long;
use strict;

my $contigFile;
my $ds = 0;
my $ups = 0;
my $usage = "usage : $0 [-num_exons] < gene locations\n";

my $num_exons = shift;

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
my @sortedGenes = ();

if(!defined($num_exons)) {
    @sortedGenes = sort {&Gene::compareTo($a, $b);} values %genes;
}
else {
    @sortedGenes = sort {&Gene::nexons_compareTo($a, $b);} values %genes;
}

foreach my $gene (@sortedGenes) {
    $gene->dump(\*STDOUT);
}
