#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use strict;

my $i;
my $id;
my $sigpinp = shift || die "$0 signalP < contigLengths\n";
open(SIGP, "<$sigpinp") || die;
my %lengths = ();
my @line = ();
while(<>) {
    @line = split(/\s+/);
    $lengths{$line[1]} = $line[0];
}

my %idsf = ();
my %idsb = ();
while(<SIGP>) {
    @line = split(/\s+/);
    if($line[0] =~ /^\s*\#/) {
	next;
    }
    if(!defined($lengths{$line[0]})) {
        next;
    }
    if(!defined($idsf{$line[0]})) {
        @{$idsf{$line[0]}} = ();
        @{$idsb{$line[0]}} = ();
        for($i = 0; $i <= $lengths{$line[0]}; $i++) {
            push(@{$idsf{$line[0]}}, 0);
            push(@{$idsb{$line[0]}}, 0);

        }
    }
    if($line[6] eq "+") {
        $idsf{$line[0]}->[$line[3]] = $line[5];
    }
    else {
        $idsb{$line[0]}->[$lengths{$line[0]} - $line[4] + 1] = $line[5];
    }
}

foreach $id (keys %idsf) {
    print ">$id\t$lengths{$id}\n";
    for($i = 1; $i <= $lengths{$id}; $i++) {
        if($idsf{$id}->[$i] != 0) {
            print "$i\t$idsf{$id}->[$i]\n";
        }
    }
    print "\/\/\n";
    for($i = 1; $i <= $lengths{$id}; $i++) {
        if($idsb{$id}->[$i] != 0) {
            print "$i\t$idsb{$id}->[$i]\n";
        }
    }

}
