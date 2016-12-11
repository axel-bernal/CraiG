#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use File::Basename;
use enum qw(STRAND_FWD STRAND_COMP);

my $match_cds = 0;
my $mismatch_tis = 0;
my $usage = "usage : $0 [--mismatch-tis] [--match-cds] annotation1 [.. annotationN]\n";

&GetOptions("--match-cds" => \$match_cds,
	    "--mismatch-tis" => \$mismatch_tis,
    );

die "Must use mismatch-tis only after match-cds is ON" if $mismatch_tis && !$match_cds;

my @genes = ();
my %clusters = ();
my $key = "";
my $value = "";
my $i = 0;
my $file;

while(scalar(@ARGV) > 0) {
    $file = shift @ARGV;
    open(FILE, "<$file") || die "Cannot open $file\n";
    push(@genes, &GenUtils::readLocsFormat(\*FILE, 0, 'TRAIN', undef));

    for my $id (keys %{$genes[-1]}) {
	my $gene = ${$genes[-1]}{$id};
	my $ngene = Gene->new(Gene => $gene);
	my $exons = $ngene->{'Exons'};

	if(scalar(@$exons) == 0) {
	    next;
	}
	
	$key = $gene->{'Phase'}."<";
	$value = [scalar(@genes) - 1, $id];

	$exons->[0]->[1] = $exons->[0]->[2] if $mismatch_tis;

	my $fexons = $ngene->{'5UTRExons'};
	my $lexons = $ngene->{'3UTRExons'};

	if(!$match_cds) {
	    for($i = 0; $i < scalar(@$lexons); $i++) {
		push(@$exons, $lexons->[$i]);
	    }
	    if(scalar(@$lexons) > 0) {
		$exons->[-1]->[2] = $exons->[-1]->[1];
	    }
	    
	    for($i = scalar(@$fexons) - 1; $i >= 0; $i--) {
		unshift(@$exons, $fexons->[$i]);
	    }
	    if(scalar(@$fexons) > 0) {
		$exons->[0]->[1] = $exons->[0]->[2];
	    }
	}

	$key .= join("_", @{$exons->[0]}[0..2]);
	for($i = 1; $i < scalar(@$exons); $i++) {
	    $key = $key.";".join("_", @{$exons->[$i]}[0..2]);
	}

#	print "key $key \n";
#	$gene->dump(\*STDOUT, undef, 1);

	my $base_file = basename($file);
	if(!defined($clusters{$key})) {
	    @{$clusters{$key}} = ();
	}
	push(@{$clusters{$key}}, $base_file);
	push(@{$clusters{$key}}, $value);
    }
    close(FILE);
}

foreach $key (keys %clusters) {
    print "//\t$key\t";
    my $c = $clusters{$key};
    for($i = 0; $i < scalar(@{$c}); $i += 2) {
	print $c->[$i];
	if($i == scalar(@{$c}) - 2) {
	    print "\n";
	}
	else {
	    print "\t";
	}
    }

    for($i = 0; $i < scalar(@{$c}); $i += 2) {
#	print "$c->[$i + 1]->[0] $c->[$i + 1]->[1]\n";
	$genes[$c->[$i + 1]->[0]]->{$c->[$i + 1]->[1]}->dump(\*STDOUT, undef, 1);
    }
}
    
