#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";

use strict;

$| = 1;

# -*- perl -*-

my $usage = "usage: $0 < bg table file > locations file\n";
my $line = "";
while(defined($line = <>) && $line =~ /(\S+)\s+(.+)$/) {
    my $id = $1;
    my @locs = split(/\s+/,$2);
    print ">$id ";
    for(my $i = 0; $i < scalar(@locs); $i +=2) {
	$locs[$i+1] -= 3 if ($i == scalar(@locs) - 2);
	print "$id","_","$locs[$i]","_",$locs[$i+1];
	if($i == scalar(@locs) - 2) {
	    print "\n";
	}
	else {
	    print  " ; ";
	}
    }
}
