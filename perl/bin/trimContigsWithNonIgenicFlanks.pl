#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $inContigs;
my $usage = "usage : $0 -ctg subContigs < igenic_locs > outContigs\n";

&GetOptions("-ctg=s" => \$inContigs,);

(defined($inContigs) && length($inContigs) > 0) || die "Contig Lengths must be specified. ".$usage;

open(CONTIGS, "<$inContigs") || die "file $inContigs not found\n";
my $c;
my $line = "";
my %lengths = ();
my $length;
my %contigs = %{&GenUtils::readLocsFormat(\*CONTIGS, 0, 'CONTIGS', undef)};
my $contig;
foreach $c (keys %contigs) {
    $contig = $contigs{$c};
    $length = $contig->{'Exons'}->[0]->[2] -
	$contig->{'Exons'}->[0]->[1] + 1;
    $lengths{$c} =$length;
}

close(CONTIGS);

my %igenics = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};

foreach my $id (keys %igenics) {
    my $ibounds = $igenics{$id}->{'Exons'}->[0];
    my $offset = 1;
    my $cid = $ibounds->[0];
    $length = $lengths{$cid};
    $contig = $contigs{$cid};

    my $cbounds = $contig->{'Exons'}->[0];
    if($cid =~/(\S+)\.([^\.])+/) {
	$cid = $1;
	$offset = $2;
    }
    
    # see if I can move the intergenic regions up or downstream
    if($ibounds->[1] != 1 && $ibounds->[1] < 150) {
	$cbounds->[1] += $ibounds->[1];
	if($offset > 1) {
	    $contig->{'Id'} = $cid.".".$cbounds->[1];
	}
    }
    if($ibounds->[2] != $length && $ibounds->[2] > $length - 150) {
	$cbounds->[1] -= $ibounds->[2];	
    }
}

foreach $c (keys %contigs) {
    $contigs{$c}->dump(\*STDOUT, undef, 1);
}
    
