#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;

my $usage = "usage : toGnuFmt.pl < scoreFile > gnu format";
my $i;
my $counter = 1;
my $line = <>;

while(defined($line) && $line =~ /^>/) {
    print $line;
    while(defined($line = <>) && $line !~ /^>/) {
	print "$counter ";
	for($i = 0; $i < 3; $i++) {
	    chomp $line; print "$line ";
	    $line = <>;
	}
	$counter += 3;
	print "\n";
    }
}
