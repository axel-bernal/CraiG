#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Genefinder;
use Contig;
use strict;

my $usage = "usage :$0 < file\n";
my @lines = ();
my $line = "";
while(defined($line = <STDIN>)) {
    chomp($line);
    my @line = split(/\s+/, $line);
    if($line =~ /^\#\s+\d+\s+\(\s*(\d+)\s*\)\s\S\s(\S+)/) {
	$line = <STDIN>; $line = <STDIN>; $line = <STDIN>;
	chomp($line);
	my $class = $line == 0 ? -1 : 1; 
	print ">$1\t$class\n$2\n";
    }
}
