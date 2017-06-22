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
my $usage = "usage : $0  < gene locations\n";

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
my %numg = ();
my $id;

foreach $id (keys %genes) {
    push(@{$numg{$genes{$id}->{'Exons'}->[0]->[0]}}, $id);
}

foreach $id (keys %numg) {
    if(scalar(@{$numg{$id}}) == 1) {
        $genes{$numg{$id}->[0]}->dump(\*STDOUT);
    }
}
