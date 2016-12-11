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
            print  ">$gene->{'Id'}.5UTR\t";
            $gene->displayLoc(\*STDOUT, $gene->{'5UTRExons'}, \%contigs);
	    print "[";
            $dist = abs($fiveU->[-1]->[2] - $exons->[0]->[1]);
            if($dist > 1) {
		print $exons->[0]->[0]."_".$exons->[0]->[1]."_".$exons->[0]->[1];
                if(scalar(@$exons) > 1) {
                    print  ">";
                }
            }
            print "\n";
        }

        if(scalar(@{$threeU}) > 0) {
            print ">$gene->{'Id'}.3UTR\t";
            $dist = abs($threeU->[0]->[1] - $exons->[-1]->[2]);
            if($dist > 1) {
                if(scalar(@$exons) > 1) {
                    print  "<";
                }
		print $exons->[-1]->[0]."_".$exons->[-1]->[2]."_".$exons->[-1]->[2];
            }
	    print "]";
            $gene->displayLoc(\*STDOUT, $gene->{'3UTRExons'}, \%contigs);
            print "\n";
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
                print STDERR $geneId, " distance last gene ", $lastGeneId, " is ", $min - $geneMax, "\n";
            }
        }
    }
}

