#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Genefinder;
use Contig;
use strict;

my $usage = "usage :$0 percentage < signal fasta file\n";
my %lines = ();
my %classes = ();
my $line = "";
my $class;
my $id;
while($line = <STDIN>) {
    chomp($line);
    if($line =~ />(\S+)\s+(\-?\d)/) {
	$id = $1; $class= $2;
    } 
    else {
	if($class == 1) {
	    print ">$id\t$class\n$line\n";
	}
	$lines{$line}->{$class} = $id;
	if(defined($lines{$line}->{1}) && defined($lines{$line}->{-1})) {
	    undef($lines{$line}->{-1});
	}
    }
}
#filtering stuff
foreach $line (keys %lines) {
    if(defined($lines{$line}->{-1})) {
	print ">$lines{$line}->{-1}\t-1\n$line\n";
    }
}

