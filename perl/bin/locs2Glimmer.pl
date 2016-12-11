#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use GenUtils;
use Gene;
use strict;
use Getopt::Long;
$| = 1;

# -*- perl -*-

my $usage = "usage: $0 < locationFile > exonsFile\n";

# loading the information
my $org = Organism->new();
my %locs = %{$org->loadGenes(\*STDIN)};
&GenUtils::writeGlimmerFormat(\%locs, \*STDOUT);
