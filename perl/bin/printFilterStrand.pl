#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $reverse = 0;
my $usage = "usage : $0 --reverse < FILTER_VALS > STRAND_FILTER_VALS\nPrints a strand of the filter. Default is forward\n";

&GetOptions("-reverse" => \$reverse,
    );

my $print_on = 1;
while(defined($_ = <STDIN>)) {    
    if($_ =~ />(\S+)/) { 
	$print_on = ($reverse ? 0 : 1);
	print $_ if !$print_on;
    } elsif($_ =~ /\/\//) { 
	$print_on = ($reverse ? 1 : 0);
    }
    print $_ if $print_on;
}
