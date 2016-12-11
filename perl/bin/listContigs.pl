#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < locations\n";

my %genes = % {&GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
foreach my $id (keys %genes) {
    print $genes{$id}->{'Exons'}->[0]->[0], "\n"; 
}
