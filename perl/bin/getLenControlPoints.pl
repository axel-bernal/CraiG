#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
use Organism;
use HistogramUtils;
use Contig;
use Getopt::Long;

my $usage = "usage : $0 contigFile exonBins intronBins igenicBins < locations\n\tGenerates length statistics\n";

my $contigFile = shift || die $usage;
my $cpExons = shift || die $usage;
my $cpIntrons = shift || die $usage;
my $cpIgenics = shift || die $usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my %genes = % { $org->loadGenes(\*STDIN)};

my %stateLens = ();
my %numCPs = ();
HistogramUtils::initGenomicRegLengths(\%stateLens, \%numCPs, 
				      $cpExons, $cpIntrons);
@{$stateLens{"INTERGENIC"}} = ();
$numCPs{"INTERGENIC"} = $cpIgenics;

for my $c (keys %contigs) {
    my $genes = $contigs{$c}->{'Features'}->{'TRAIN'};

    foreach my $id (keys %$genes) {
	HistogramUtils::extractGenomicRegLengths(\%stateLens, 
						 $genes->{$id}->{'Exons'},
						 $genes->{$id}->{'Strand'});
    }
    
    # dealing with igenics
    my $array = &GenUtils::markIGenics($contigs{$c}, 0, 1, 'TRAIN');
    my $interBeg = -1;
    my $interEnd = -1;

    for(my $i = 1; $i <= $contigs{$c}->{'Len'}; $i++) {
        if($array->[$i] == 0 && $interBeg == -1) {
            $interBeg = $i;
        }
        if($i > 1) {
            if(($array->[$i] == 1 && $array->[$i - 1] == 0) || $i == $contigs{$c}->{'Len'}) {
                $interEnd = $i - 1;
            }
        }

        if($interBeg != -1 && $interEnd != -1) {
            push(@{$stateLens{"INTERGENIC"}}, [$c, $interBeg, $interEnd]);
            $interBeg = -1;
            $interEnd = -1;
        }
    }
}

foreach my $keyw (keys %stateLens) {
#    print STDERR $keyw, "\n";
    my $st_lens = $stateLens{$keyw};
    my $num_cps = $numCPs{$keyw};
    my $cps = HistogramUtils::getLengthControlPoints($st_lens, $num_cps);
    HistogramUtils::printLengthControlPoints(\*STDOUT, $st_lens, $cps, $keyw);
}
