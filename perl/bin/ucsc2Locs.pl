#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use GenUtils;
use Gene;
use strict;
use Getopt::Long;
$| = 1;

# -*- perl -*-

my $usage = "usage: $0 -offset offset < gff file > locations file\n";
my $offset = 0;
&GetOptions("--offset=i" => \$offset);
$offset >= 0 || die "offset must be >= 0. $usage";
# loading the information
my %genes = % { &GenUtils::readUCSCFormat(\*STDIN, 1, $offset)};
&GenUtils::displayGeneList(\%genes, \*STDOUT);
