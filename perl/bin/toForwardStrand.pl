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
my $locs;

my $usage = "usage : $0 -ctg contigFile -f locations. All files in fasta format\n";

&GetOptions("-f=s" => \$locs,
	    "-ctg=s" => \$contigFile,
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(FASTA, "<$locs") || die "$locs\n";
open(MOD_CONTIG, ">$contigFile.fwd") || die "$contigFile.fwd\n";
open(MOD_FASTA, ">$locs.fwd") || die "$locs.fwd\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*FASTA)};
foreach my $cid (keys %contigs) {
    my $c = $contigs{$cid};
    my $geneOnComp = 0;
    my $geneOnFwd = 0;
    foreach my $gid (keys % {$c->{'Features'}->{'TRAIN'}}) {
	my $g = $c->{'Features'}->{'TRAIN'}->{$gid};
	foreach my $exon (@{$g->{'Exons'}}) {
	    $exon->[0] = "$c->{'Id'}.$g->{'Strand'}";
	}
	if($g->{'Strand'} == STRAND_COMP) {
	    @{$g->{'Exons'}} = reverse(@{$g->{'Exons'}});
	    foreach my $exon (@{$g->{'Exons'}}) {
		my $tmp = $exon->[1];
		$exon->[1] = $exon->[2];
		$exon->[2] = $tmp;
	    }
	    $geneOnComp = 1;
	}
	else {
	    $geneOnFwd = 1;
	}
	$g->displayHeader(\*MOD_FASTA);
	print MOD_FASTA "\n";
	if($geneOnFwd) {
	    print MOD_CONTIG ">$c->{'Id'}.".STRAND_FWD."\n$c->{'Seq'}\n";
	}
	if($geneOnComp) {
	    print MOD_CONTIG ">$c->{'Id'}.".STRAND_COMP."\n".GenUtils::revComp($c->{'Seq'})."\n";
	}
    }
}
