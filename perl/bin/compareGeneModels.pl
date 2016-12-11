#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $org = Organism->new();

my $usage = "usage : $0 gene_file_1 gene_file_2 contigs > differences\n";

my $file_a = shift || die "must provide gene_file_1 as parameter\n.".$usage;
my $file_b = shift || die "must provide gene_file_2 as parameter\n.".$usage;
my $contig_file = shift || die "must provide contigs\n".$usage;

open(CONTIGS, "<$contig_file") || die "Can't open contig_file.\n".$usage;

my %contigs = % {$org->loadContigs(\*CONTIGS)};

open(GENES_A, "<$file_a") || die "Can't open gene_file_1.\n".$usage;
open(GENES_B, "<$file_b") || die "Can't open gene_file_2.\n".$usage;

my %genes_a = % { $org->loadGenes(\*GENES_A, 'GENES_A')};
my %genes_b = % { $org->loadGenes(\*GENES_B, 'GENES_B')};
my @genes_a = sort {&Gene::compareTo($a, $b);} (values %genes_a);
my @genes_b = sort {&Gene::compareTo($a, $b);} (values %genes_b);

my @maps = (0,0,0);
my $toth = 0;
my $totl = 0;
my $num_genes = scalar(@genes_a);
for(my $i = 0; $i < $num_genes; $i++) {
    my $gene_a = $genes_a[$i];
    my $gene_b = $genes_b[$i];

    my $lowest_a = $gene_a->lowestCoord();
    my $highest_a = $gene_a->highestCoord();
    my $lowest_b = $gene_b->lowestCoord();
    my $highest_b = $gene_b->highestCoord();
    my $diffl = abs($lowest_a - $lowest_b);
    my $diffh = abs($highest_a - $highest_b);
    my $cid = $gene_a->{'Contig'}->{'Id'};
    $cid =~ /(\S+)\.(\d+)$/;
    print $diffh + $diffl, " $diffl $diffh (", $lowest_a, ", ", $highest_a, ") (",$lowest_a+$2-1,", ",$highest_a+$2-1,") (",$lowest_b+$2-1,",", $highest_b+$2-1,") $1:$2..", $2 + $gene_a->{'Contig'}->{'Len'},"\n";
    if(!$diffl) {
	$maps[0]++;
    }
    if(!$diffh) {
	$maps[1]++;
	if(!$diffl) {
	    $maps[2]++;
	}
    }
    $totl += $diffl;
    $toth += $diffh;
}

print "5prime matches ", $maps[0]/$num_genes, "\n";
print "3prime matches ", $maps[1]/$num_genes, "\n";
print "both matches ", $maps[2]/$num_genes, "\n";
print "avg difference 5prime end ", $totl/$num_genes, "\n";
print "avg difference 3prime end ", $toth/$num_genes, "\n";
