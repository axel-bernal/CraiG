#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Genefinder;
use Contig;
use strict;

my $usage = "usage :$0 percentage < file\n";
my $perc = shift || die "$usage";
my @lines = ();
my $line = "";
while($line = <STDIN>) {
    chomp($line);
    push(@lines, $line);
}

my $num = 0;
my %seen = ();
while($num < (scalar(@lines)*$perc/100.0)) {
    my $r = int (rand (scalar(@lines) - 1));
    while(defined($seen{$r})) {
	$r = int (rand (scalar(@lines) - 1));	
    }
    print "$lines[$r]\n";
    $seen{$r} = 1;
    $num++;
}
