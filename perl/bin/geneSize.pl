#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $usage = "usage : $0 < locations > sizes\n";

my $line;
while(defined($line = <STDIN>) && ($line =~ /^>(\S+)\s+(.+)\s+?(\d+)?b?p?$/)) {
    print "$1";
    my (undef, undef, undef, $locs, undef, undef) = @{&GenUtils::readContigLocation($2)};
    my $size = 0;
    foreach my $loc (@$locs) {
        $size += abs($loc->[2] - $loc->[1]) + 1;
    }
    my $threeP = ($size % 3 == 0) ? 1 : 0;
    print " $size ($threeP)\n";
}
