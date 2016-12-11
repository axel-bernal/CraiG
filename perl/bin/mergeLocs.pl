#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < locations > mergedLocs\n";

my $line;
while(defined($line = <STDIN>) && ($line =~ /^>(\S+)\s+(.+)\s+?(\d+)?b?p?$/)) {
    print ">$1\t";
    my @locs = @{&GenUtils::readLocExons($2)};
    my $lexon = shift @locs;

    foreach my $loc (@locs) {
        if($lexon->[0] eq $loc->[0] && $loc->[1] == $lexon->[2] + 1) {
            $lexon->[2] = $loc->[2];
        }
        else {
            print "$lexon->[0]_$lexon->[1]_$lexon->[2];";
            $lexon->[0] = $loc->[0];
            $lexon->[1] = $loc->[1];
            $lexon->[2] = $loc->[2];
        }
    }
    print "$lexon->[0]_$lexon->[1]_$lexon->[2]\n";
}
