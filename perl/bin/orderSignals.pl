#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Genefinder;
use Contig;
use strict;

my $usage = "usage :$0 percentage < signal fasta file\n";
my %strands = ();
my %lines = ();
my $line;
while(defined($line = <STDIN>) && $line =~ />(\S+)_(\d)_(\d+)\s+(\-?\d)/) {
    $strands{$2}->{$1}->{$3} = $4; 
    $line = <STDIN>;
    $lines{$2}->{$1}->{$3} = $line;
} 

#filtering stuff
foreach my $strand (keys %strands) {
    foreach my $contig (sort {$a cmp $b} keys %{$strands{$strand}}) {
	foreach my $pos (sort {$a <=> $b} keys %{$strands{$strand}->{$contig}}) {
	    print ">$contig"."_"."$strand"."_"."$pos\t", $strands{$strand}->{$contig}->{$pos}, "\n";
	    print $lines{$strand}->{$contig}->{$pos};
	}
    }
}
