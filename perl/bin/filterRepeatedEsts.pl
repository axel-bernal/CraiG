#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $contigFile;
my $est;
my $last_est;
my $usage = "usage : $0 -ctg contigFile < locations > filter\n";

&GetOptions(
            "-ctg=s" => \$contigFile,
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my %ests = % { $org->loadGenes(\*STDIN)};

my @ests = sort {&Gene::compareTo($a, $b);} (values %ests);
if(scalar(@ests) == 0) {
    exit;
}

my $j = 0;
my $length;
my $last_length;

for(my $i = 1; $i < scalar(@ests); $i++) {
    $last_est = $ests[$j];
    $est = $ests[$i];
    if(equalTo($est, $last_est)) { 
        if(scalar(@{$last_est->{'Exons'}}) == 1) { #overlapping single-exons
            if($est->{'Strand'} == STRAND_FWD) {
                $est->{'Exons'}->[0]->[1] = Utils::min($last_est->{'Exons'}->[0]->[1], $est->{'Exons'}->[0]->[1]);
                $est->{'Exons'}->[0]->[2] = Utils::max($last_est->{'Exons'}->[0]->[2], $est->{'Exons'}->[0]->[2]);
            }
            else {
                $est->{'Exons'}->[0]->[1] = Utils::max($last_est->{'Exons'}->[0]->[1], $est->{'Exons'}->[0]->[1]);
                $est->{'Exons'}->[0]->[2] = Utils::min($last_est->{'Exons'}->[0]->[2], $est->{'Exons'}->[0]->[2]);
            }
            splice(@ests, $j, 1);
        }
        else {
            $length = $est->highestCoord() - $est->lowestCoord();
            $last_length = $last_est->highestCoord() - $last_est->lowestCoord();
            if($last_length > $length) {
                splice(@ests, $i, 1);
            }
            else {
                splice(@ests, $j, 1);
            }
        }
        $i--;
    }
    else {
        $j = $i;
    }
}

foreach $est (@ests) {
    if(goodAlignment($est) && scalar(@{$est->{'Exons'}}) > 0) {
        $est->dump(\*STDOUT, undef, 1);
    }
}


sub equalTo {
    my ($est1, $est2) = @_;

    if(scalar(@{$est1->{'Exons'}}) != scalar(@{$est2->{'Exons'}})) {
        return 0;
    }
    for(my $i = 0; $i < scalar(@{$est1->{'Exons'}}); $i++) {
        my $exon1 = $est1->{'Exons'}->[$i];
        my $exon2 = $est2->{'Exons'}->[$i];
        if($exon1->[0] ne $exon2->[0]) {
            return 0;
        }
        if($i != 0 && $exon1->[1] != $exon2->[1]) {
            return 0;
        }
        if($i != scalar(@{$est1->{'Exons'}}) - 1 && $exon1->[2] != $exon2->[2]) {
            return 0;
        }
    }
    if(scalar(@{$est1->{'Exons'}}) == 1) {
        if(($est1->lowestCoord() >= $est2->lowestCoord() && $est1->lowestCoord() <= $est2->highestCoord()) ||
           ($est1->highestCoord() >= $est2->lowestCoord() && $est1->highestCoord() <= $est2->highestCoord())
           ) {
            return 1;
        }
        else {
            return 0;
        }
    }
    return 1;
}

sub goodAlignment  {
    my ($est) = @_;
    my $donor;
    my $acceptor;
    my $strand = $est->{'Strand'};
    my $e_i;

    for(my $i = 0; $i < scalar(@{$est->{'Exons'}}); $i++) {
        $e_i = $est->{'Exons'}->[$i];
        if($i > 0) {
            my $e_i_1 = $est->{'Exons'}->[$i - 1];
            $donor = substr($est->{'Contig'}->getSequence(0, $e_i_1->[1], $e_i_1->[2], 2, 1, $strand), abs($e_i_1->[2] - $e_i_1->[1]) + 1, 2);
            $acceptor = substr($est->{'Contig'}->getSequence(2, $e_i->[1], $e_i->[2], 0, 1, $strand), 0, 2);
            
            if($donor !~ /gt/i) {
                return 0;
            }
            if($acceptor !~ /ag/i) {
                return 0;
            }
        }
    }
    if(scalar(@{$est->{'Exons'}}) == 1) {
        return 0;
    }
    return 1;
}
