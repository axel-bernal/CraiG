#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Gene;
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $ds = 0;
my $ups = 0;
my $usage = "usage : $0 -ctg contigFile < locations > revcompLocations\n";

&GetOptions("-ups=i" => \$ups,
	    "-ds=i" => \$ds,
	    "-ctg=s" => \$contigFile,
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();
my $gene;
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};
foreach my $id (keys %genes) {
    $gene = $genes{$id};
    foreach my $exon (@{$gene->{'Exons'}}) {
	$exon->[2] = $gene->{'Contig'}->{'Len'} - $exon->[2] + 1;
	$exon->[1] = $gene->{'Contig'}->{'Len'} - $exon->[1] + 1;
    }
}

my @genes = values %genes;
my @sortedGenes = sort {&Gene::compareTo($a, $b);} (@genes);
#       print  scalar(@genes),"\n";
foreach $gene (@sortedGenes) {
    $gene->displayHeader(\*STDOUT, \%contigs);
    print "\n";
}
