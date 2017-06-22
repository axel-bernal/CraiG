#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 [-l] -ctg subcontig_info < score_filter > new_score_filter\n";

my $subcontig = "";
my $subContigsAreLocs = 0;

&GetOptions("-l" => \$subContigsAreLocs,
            "-ctg=s" => \$subcontig,
            );

open(SUBCONTIG, "$subcontig") ||  die $usage;

my $scLengths;
my $contigRanges;

if($subContigsAreLocs) {
    my $line;
    ($contigRanges, $scLengths) = @{&GenUtils::loadSubContigLocs(\*SUBCONTIG)};
}
else {
    ($contigRanges, $scLengths) = @{&GenUtils::loadSubContigSeqs(\*SUBCONTIG)};    
}

close(SUBCONTIG);

my $line =  <STDIN>;

while(defined($line) && $line  =~ />(\S+)\s(\d+)/) {
    my $cid = $1;
    my $clen = $2;
    my $strand = STRAND_FWD;
    
    if(!defined($contigRanges->{$cid})) {
        while(defined($line = <STDIN>) && $line  !~ />\S+\s\d+/) { }
        next;
    }
    
    my %subContigScores = ();    
    my $index = 0;
    my $range = $contigRanges->{$cid}->[$index];
    my $length =  $range->[3] - $range->[2] + 1;
    my $scores = GenUtils::createFilterValArrays(\%subContigScores, $range->[0], $length);
#    print "$cid $clen $length $range->[0] $index\n";    

    while(defined($line = <STDIN>) && $line  !~ />(\S+)\s(\d+)/) {

        my $sigPos1;
	my $sigPos2;
        my $sigScore;
        if($line =~ /^([\d,\.]+)\s+([\-,\d,.]+)\s*$/) {
            $sigPos1 = $sigPos2 = $1;
            $sigScore = $2;
	    if($1 =~ /^([^\.]+)\.\.(\d+)$/) {
		$sigPos1 = $1;
		$sigPos2 = $2;
	    }
	    for(my $sigPos = $sigPos1; $sigPos <= $sigPos2; $sigPos++) {
		if($strand == STRAND_FWD) {
		    while($sigPos > $range->[3] && 
			  $index < scalar(@{$contigRanges->{$cid}}) - 1) {
			$index++;
			$range = $contigRanges->{$cid}->[$index];
			$length =  $range->[3] - $range->[2] + 1;
			$scores = GenUtils::createFilterValArrays(\%subContigScores, $range->[0], $length);
#		    print "creating $range->[0]\n";
		    }
#                print "$cid $strand $sigPos $index ", $range->[2], " ", $range->[3], " $sigScore\n";
		    
		    if($sigPos < $range->[2]) {
			
		    }
		    elsif($sigPos <= $range->[3]) {
#		    print "inserted $range->[0]\n";
			$scores->[$strand]->[$sigPos - $range->[2] + 1] = $sigScore;                
		    }
		}
		else {
		    while($sigPos > $clen - $range->[2] + 1 &&
			  $index > 0) {
			$index--;
			$range = $contigRanges->{$cid}->[$index];
			$length =  $range->[3] - $range->[2] + 1;
			$scores = $subContigScores{$range->[0]};
#		    $scores = GenUtils::createFilterValArrays(\%subContigScores, $range->[0], $length);
		    }
		    
#                print "$cid $strand $sigPos $index ", $clen - $range->[3] + 1, " ", $clen - $range->[2] + 1, " $sigScore\n";
		    if($sigPos < $clen - $range->[3] + 1) {
			
		    }
		    elsif($sigPos <= $clen - $range->[2] + 1) {
#                    print "inserted $range->[0]\n";
			$scores->[$strand]->[$sigPos - $clen + $range->[3]] = $sigScore;
		    }
		}
	    }
        }
        elsif($line =~ /\/\//) {
            $strand = STRAND_COMP;
        }
    }

    while($index < scalar(@{$contigRanges->{$cid}}) - 1) {
	$index++;
	$range = $contigRanges->{$cid}->[$index];
	$length =  $range->[3] - $range->[2] + 1;
	$scores = GenUtils::createFilterValArrays(\%subContigScores, $range->[0], $length);
    }
    
#    print "contig $cid ", scalar keys %subContigScores,"\n";;    
    GenUtils::printScoreFilterVals(\%subContigScores, \*STDOUT, $strand);
    foreach $subcontig (keys %subContigScores) {
	undef $subContigScores{$subcontig};
    }
}
