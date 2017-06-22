#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 < signal_filter > aggregated_signal_filter\n";

my $line = <STDIN>;
my %filterVals = ();

while(1) {
    my $cid = GenUtils::loadFilterVals4Id(\*STDIN, \$line, \%filterVals, 0);
    if($cid eq "") {
	last;
    }
}

GenUtils::printScoreFilterVals(\%filterVals, \*STDOUT, STRAND_COMP);
print STDERR "Done\n";
