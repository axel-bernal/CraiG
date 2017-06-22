#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use GenUtils;
use Getopt::Long;

my $usage = "usage : $0 < genes\n";

my %genes = % {&GenUtils::readWeightedFastaFormat(\*STDIN)};
my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);

my $olapStrand = [0,0,0];
my %olapLengths = ();

if(!scalar(@genes)) {
    exit;
}

my $length;
my @sgene = ($genes[0]);
my $sgene_end = $genes[0]->highestCoord();

for(my $i = 1; $i < scalar(@genes); $i++) {
    my $lgene = $genes[$i - 1];
    my $rgene = $genes[$i];

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
	my $maxId = pickLongestTranscript(\@sgene);
	if(defined($maxId)) {
	    $sgene[$maxId]->displayHeader(\*STDOUT, undef, 1);
	    print "\t$sgene[$maxId]->{'Exons'}->[0]->[3]";
	    for(my $j = 1; $j < scalar(@{$sgene[$maxId]->{'Exons'}}); $j++) {
		print ";", $sgene[$maxId]->{'Exons'}->[$j]->[3];
	    }
	    print "\n";
	}
	
	@sgene = ($rgene);
	$sgene_end = $rgene->highestCoord();
    }
}

if($#sgene >= 0) {
    my $maxId = pickLongestTranscript(\@sgene);
    if(defined($maxId)) {
	$sgene[$maxId]->displayHeader(\*STDOUT, undef, 1);
	print "\t$sgene[$maxId]->{'Exons'}->[0]->[3]";
	for(my $j = 1; $j < scalar(@{$sgene[$maxId]->{'Exons'}}); $j++) {
	    print ";", $sgene[$maxId]->{'Exons'}->[$j]->[3];
	}
	print "\n";
    }
}	

sub pickLongestTranscript {
    my($sgene) = @_;
    my $maxCodingLength = 0;
    my $maxLength = 0;
    my $maxId = undef;	    
    
    for(my $j = 0; $j < scalar(@$sgene); $j++) {
	my $transcript = $sgene[$j];
	my $exons = $transcript->{'Exons'};
	if(scalar(@$exons) == 0) {
	    next;
	}
	
	my $transcriptLength = $transcript->highestCoord() - 
	    $transcript->lowestCoord() + 1;
	
	if($transcriptLength > $maxLength) {
	    $maxId = $j;
	    $maxLength = $transcriptLength;
	}
	elsif($transcriptLength == $maxLength) {
	    my $codingLength = 0;
	    for(my $j = 0; $j < scalar(@$exons); $j++) {
		$codingLength += abs($exons->[$j]->[2] - $exons->[$j]->[1]) + 1;
	    }
	    if($codingLength > $maxCodingLength) {
		$maxId = $j;
		$maxCodingLength = $codingLength;
	    }
	}
    }
    
    return $maxId;
}
