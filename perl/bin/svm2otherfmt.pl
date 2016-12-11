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
while($line = <STDIN>) {
    chomp($line);
    my @line = split(/\s+/, $line);
    my $class = shift @line;
    foreach my $elem (@line) {
	my($ind, $val) = split(/\:/,$elem);
	print "$val ";
    }
    if($class == -1){
	print "0\n";
    }
    else {
	print "1\n";
    }
}
