#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
use Organism;
use Contig;

my $usage = "usage : $0 numBins < sorted_list\n\tGenerates statistics for the list(distribution)\n";

my $numBins = shift || die $usage;
my $line;
my @list = ();

while(defined($line = <STDIN>)) {
    chomp $line;
    push(@list, $line);
}
printLengthControlPoints(\@list, $numBins);

sub printLengthControlPoints {
    my($array, $cps) = @_;
    my $offset = 0;
    my $e;

    for(my $i = 0; $i < $cps; $i++) {
        $e = $array->[int($offset)];
        print "[", int($offset),"] -> ", $e,"\n";
        $offset += 1.0*scalar(@$array)/$cps;
    }
    $e =  $array->[-1];
    print "[", scalar(@$array),"] -> ", $e, "\n";
}
