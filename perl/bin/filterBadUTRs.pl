#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);


my $usage = "usage : $0 contig < locations\n";

my $org = Organism->new();
my $gene;
my $contigFile = shift;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my %genes = % { $org->loadGenes(\*STDIN)};

my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);
if(scalar(@genes) == 0) {
    exit;
}

my $dist;
my $strand;

for my $c (keys %contigs) {
    my @array = ();
    my $i;
    my $gene;

    my @genes = sort {&Gene::compareTo($a, $b);} (values %{ $contigs{$c}->{'Features'}->{'TRAIN'}});
#    print "$c ", $contigs{$c}->{'Len'}," ", scalar(@genes), "\n";

    my $lastGeneId = "";
    my $geneMax = 1;
    my $geneId = "";
    
    for $gene (@genes) {
        $strand = $gene->{'Strand'};
        my $min = 0;
        my $max = 0;

        my $exons = $gene->{'Exons'};
        my $fiveU = $gene->{'5UTRExons'};
        my $threeU = $gene->{'3UTRExons'};
        if(scalar(@{$fiveU}) > 0) {
            $dist = abs($fiveU->[-1]->[2] - $exons->[0]->[1]);
            if($dist > 1) {
                print $gene->{'Id'}, " distance 5UTR = ", $dist, "\n";
            }
        }

        if(scalar(@{$threeU}) > 0) {
            $dist = abs($threeU->[0]->[1] - $exons->[-1]->[2]);
            if($dist > 1) {
                print $gene->{'Id'}, " distance 3UTR = ", $dist, "\n";
            }
        }

        if($gene->{'Id'} =~ /(\S+)\|([^\-]+)$/) {
            $geneId = $1;
        }
        else {
            $geneId = $gene->{'Id'};
        }
        
        $max = $gene->highestCoord();
        
        if($lastGeneId ne "" && $geneId eq $lastGeneId) { # overlap allowed between splice variants
            $geneMax = Utils::max($geneMax, $max);
        }
        else {
            $min = $gene->lowestCoord();
            if($min < $geneMax) {
                print $geneId, " distance last gene ", $lastGeneId, " is ", $min - $geneMax, "\n";
            }
        }
    }
}

