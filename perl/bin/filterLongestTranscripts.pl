#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use GenUtils;
use Getopt::Long;

my $usage = "usage : $0 < genes\n";

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
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
	    $sgene[$maxId]->dump(\*STDOUT, undef,  1);
	}
	
	@sgene = ($rgene);
	$sgene_end = $rgene->highestCoord();
    }
}

if($#sgene >= 0) {
    my $maxId = pickLongestTranscript(\@sgene);
    if(defined($maxId)) {
	$sgene[$maxId]->dump(\*STDOUT, undef,  1);
    }
}	

sub pickLongestTranscript {
    my($sgene) = @_;
    my $maxId = undef;	    
    my $maxCodingLength = 0;
    my $maxLength = 0;
    
    for(my $j = 0; $j < scalar(@$sgene); $j++) {
	my $transcript = $sgene[$j];
#	print STDERR "examining $transcript->{'Id'} $j\n";
	my $exons = $transcript->{'Exons'};
	if(scalar(@$exons) == 0) {
	    next;
	}
	
	my $transcriptLength = $transcript->highestCoord() - 
	    $transcript->lowestCoord() + 1;
	
	if($transcriptLength > $maxLength) {
	    $maxId = $j;
	    $maxLength = $transcriptLength;
	    my $codingLength = 0;
	    for(my $k = 0; $k < scalar(@$exons); $k++) {
		$codingLength += abs($exons->[$k]->[2] - $exons->[$k]->[1]) + 1;
	    }
	    $maxCodingLength = $codingLength;
	}
	elsif($transcriptLength == $maxLength) {
	    my $codingLength = 0;
	    for(my $k = 0; $k < scalar(@$exons); $k++) {
		$codingLength += abs($exons->[$k]->[2] - $exons->[$k]->[1]) + 1;
	    }
#	    print STDERR "equal lengths $transcript->{'Id'} $maxCodingLength $codingLength\n";
	    if($codingLength > $maxCodingLength) {
		$maxId = $j;
		$maxCodingLength = $codingLength;
	    }
	}
    }
#    print STDERR "choosing $maxId\n";
    return $maxId;
}
