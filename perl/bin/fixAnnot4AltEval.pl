#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use GenUtils;
use Getopt::Long;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 < genes\n";

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);

my @str_genes = ([], []);

if(!scalar(@genes)) {
    exit;
}

my $id = 1;

for(my $i = 0; $i < scalar(@genes); $i++) {
    push(@{$str_genes[$genes[$i]->{'Strand'}]}, $genes[$i]);
}

for(my $s = 0; $s < 2; $s++) {
    my $genes = $str_genes[$s];
    my $length;
    my @sgene = ($genes->[0]);
    my $sgene_end = $genes->[0]->highestCoord();
    
    for(my $i = 1; $i < scalar(@$genes); $i++) {
	my $lgene = $genes->[$i - 1];
	my $rgene = $genes->[$i];
	
	my $b = $rgene->compareTo($lgene);
	
	my $olap = 0;
	if($lgene->contigId() eq $rgene->contigId()) {
	    $olap = $sgene_end - $rgene->lowestCoord();
	}
	
	
	if($olap > 0) {
#	print "overlap $sgene_end ", $rgene->lowestCoord(), " $lgene->{'Id'} $rgene->{'Id'}\n";
	    push(@sgene, $rgene);
	    $sgene_end = ($sgene_end > $rgene->highestCoord()) ?
		$sgene_end : $rgene->highestCoord();
	}
	else {
	    dumpSgene($id++, \@sgene);
	    @sgene = ($rgene);
	    $sgene_end = $rgene->highestCoord();
	}
    }
    
    dumpSgene($id++, \@sgene);
}

sub dumpSgene {
    my($id, $sgene) = @_;

    for(my $j = 0; $j < scalar(@$sgene); $j++) {
	my $gene = $sgene->[$j];
#	print STDERR "examining $transcript->{'Id'} $j\n";
	my $exons = $gene->{'Exons'};
	if(scalar(@$exons) == 0) {
	    next;
	}

	print  ">gene_$id|$gene->{'tId'}\t";
	$gene->displayLoc(\*STDOUT, $gene->{'Exons'});
	print "\n";
    }
}
