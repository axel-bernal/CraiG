#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Gene;
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $contigFile;
my $strand = 0;
my $usage = "usage : $0 -strand <strand> -ctg <contigFile> < locationsFile > revcompLocations\n";

&GetOptions("-strand=i" => \$strand,
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
    if($gene->{'Strand'} != $strand) {
        foreach my $exon (@{$gene->{'Exons'}}) {
            $exon->[2] = $gene->{'Contig'}->{'Len'} - $exon->[2] + 1;
            $exon->[1] = $gene->{'Contig'}->{'Len'} - $exon->[1] + 1;
        }
    }
    $gene->displayHeader(\*STDOUT, \%contigs);
    if($gene->{'Strand'} != $strand) {
        print "\n", GenUtils::revComp($gene->{'Contig'}->{'Seq'}), "\n";
    }
    else {
        print "\n", $gene->{'Contig'}->{'Seq'}, "\n";
    }
}
