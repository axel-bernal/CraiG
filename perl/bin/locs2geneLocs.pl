#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < locations\n";
my $line;
my $id;

while(defined($line = <STDIN>) && ($line =~ /^>(\S+)\s+(.+)\s+?(\d+)?b?p?$/)) {
    my $id = $1;
    my @locs = @{&GenUtils::readLocExons($2)};
    if(!scalar(@locs)) {
        next;
    }

    my $contig = "";
    my $part = 1;

    foreach my $loc (@locs) {
        if($loc->[0] ne $contig) {
            if($contig ne "") {
                print "\n";
            }
            print ">$id\@part:$part\t";
            $part++;
        }
        else {
            if($contig ne "") {
                print ";";
            }
        }
        print join("_", @$loc);
        $contig = $loc->[0];
    }
    print "\n";
}
