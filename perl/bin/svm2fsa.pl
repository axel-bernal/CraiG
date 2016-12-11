#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Genefinder;
use Contig;
use strict;

my $usage = "usage :$0 <  svmFile > fastaFile\n";
my $line = "";

my(%encoding);

$encoding{'0000'} = 'X';
$encoding{'0001'} = 'A';
$encoding{'0010'} = 'C';
$encoding{'0100'} = 'T';
$encoding{'1000'} = 'G';

my @linesSVM = ();
while($line = <STDIN>) {
    chomp($line);
    my @line = split(/\s+\d+\:/, $line);
    push(@linesSVM, \@line);
}

for(my $i = 0; $i < scalar(@linesSVM); $i++) {
    my $class = shift @{$linesSVM[$i]};
    print ">id_$i\t$class\n";
#    print join(" ",@{$linesSVM[$i]}), "\n";
    for(my $j = 0; $j <= scalar(@{$linesSVM[$i]}) - 4; $j += 4) {
	print $encoding{$linesSVM[$i]->[$j].$linesSVM[$i]->[$j+1].$linesSVM[$i]->[$j+2].$linesSVM[$i]->[$j+3]};
    }
    print "\n";
}
