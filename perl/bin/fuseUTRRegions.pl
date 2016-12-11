#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use GenUtils;
use enum qw(STRAND_FWD STRAND_COMP);
my %genes = %{&GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
&GenUtils::fuseUtrRegions(\%genes);
&GenUtils::displayGeneList(\%genes, \*STDOUT);
