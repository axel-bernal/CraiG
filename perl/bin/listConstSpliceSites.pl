#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $usage = "usage : $0 < locations\n";
my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
my %donors;
my %acceptors;

foreach my $id (keys %genes) {
    my $strand = $genes{$id}->{'Strand'};
    my $contig = $genes{$id}->{'Exons'}->[0]->[0];

    for(my $i = 0; $i < scalar(@{$genes{$id}->{'Exons'}}); $i++) {
        if(scalar(@{$genes{$id}->{'Exons'}}) == 1) {
            next;
        }

        my $acceptorP = $genes{$id}->{'Exons'}->[$i]->[1] - 2;
        my $donorP = $genes{$id}->{'Exons'}->[$i]->[2] + 1;

        if($strand == STRAND_COMP) {
            $acceptorP = $genes{$id}->{'Exons'}->[$i]->[1] + 2;
            $donorP = $genes{$id}->{'Exons'}->[$i]->[2] - 1;
        }
        
        if($i != 0) {
            ${$acceptors{$contig}{$strand}}{$acceptorP}++;
        }

        if($i != scalar(@{$genes{$id}->{'Exons'}}) - 1) {
            ${$donors{$contig}{$strand}}{$donorP}++;
        }
    }
}

foreach my $c (keys %donors) {
    foreach my $strand (keys %{$donors{$c}}) {
        foreach my $p (keys %{$donors{$c}{$strand}}) {
            if(${$donors{$c}{$strand}}{$p} > 1) {
                print "$c 2 $strand $p\n";
            }
        }
    }
}

foreach my $c (keys %acceptors) {
    foreach my $strand (keys %{$acceptors{$c}}) {
        foreach my $p (keys %{$acceptors{$c}{$strand}}) {
            if(${$acceptors{$c}{$strand}}{$p} > 1) {
                print "$c 3 $strand $p\n";
            }
        }
    }
}
