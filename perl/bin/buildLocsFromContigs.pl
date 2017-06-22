#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

my $usage = "$0 seqLens_in_asc_order contigLens_in_asc_order > locs\n";
my $seqLens = shift || die $usage;
my $cLens = shift || die $usage;

open(SEQL, "<$seqLens") || die "Can't open $seqLens\n";
open(CONTIGL, "<$cLens") || die "Can't open $cLens\n";

my @cLens;

my $line;
while(defined($line = <CONTIGL>) && $line =~ /^(\d+)\s+(\S+)$/) {
    push(@cLens, [$1, $2]);
}

my $currCtg = 0;
my $currEndP = 0;

while(defined($line = <SEQL>) && $line =~ /^(\d+)\s+(\S+)$/) {
    my $seqLen = $1;
    my $seqId = $2;
    if($currCtg >= scalar(@cLens))  {
        last;
    }
    print ">$seqId\t";
    while($seqLen && $currCtg < scalar(@cLens)) {
        my $remnant =$cLens[$currCtg]->[0] - $currEndP;
        if(!$remnant) {
            $currCtg++;
            $currEndP = 0;
            next;
        }

        if($remnant < $seqLen) {
            $seqLen -= $remnant;
            print $cLens[$currCtg]->[1],"_", $currEndP + 1, "_", $cLens[$currCtg]->[0], ";";
            $currEndP = 0;
            $currCtg++;

        }
        else {
            print $cLens[$currCtg]->[1],"_", $currEndP + 1, "_", $seqLen + $currEndP, "\n";
            $currEndP += $seqLen;
            $seqLen = 0;
        }
    }
}
      
