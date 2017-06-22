#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Genefinder;
use Contig;
use strict;

my $usage = "usage :$0 < GenSplicer_file\n";


my $line = "";
my $begin = 1;
my $oldId = "";
my @line = ();
while(defined($line  = <STDIN>)) { 
    chomp($line);
    if($line eq "") {
	next;
    }
    $line =~ /(\S+)\s+/;
    if(scalar(@line)) {
	if($oldId ne "" && $oldId ne $1) {
	    $line[2] -= 3;
	}
	print join("_", @line);
    }
    @line = split(/\s+/, $line);
    if($line[0] ne $oldId) {
	if($oldId ne "") {
	    print "\n";
	}
	$oldId = $line[0];
	$begin = 1;
    }
    else {
	print " ; ";
    }
    if($begin == 1) {
	print ">$line[0]"."_"."$begin ";
    }
    $begin++;
}
