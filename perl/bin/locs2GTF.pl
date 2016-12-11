#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use GenUtils;
use Contig;
use Gene;
use strict;
use Getopt::Long;
$| = 1;

# -*- perl -*-
my $locs = "";
my $contigs = "";

my $usage = "usage: $0 contig_file < locationFile > exonsFile\n";
$contigs = shift || die "$usage";
open(CONTIGS, "<$contigs") || die "Can't open file $contigs\n";
# loading the information
#my %contigs = %{$org->loadContigs(\*CONTIGS)};
my %locs = %{&GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
&GenUtils::writeGTFFormat(\%locs, \*STDOUT);
