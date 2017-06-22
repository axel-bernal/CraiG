#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Gene;
use GenUtils;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $contigFile;
my $ds = 0;
my $ups = 0;
my $usage = "usage : $0 annotations start_score\n";

my $annotFile = shift || die $usage;
my $scoreFile = shift || die $usage;

open(ANNOT, "$annotFile") || die $usage;
open(SCORE, "$scoreFile") || die $usage;

my %genes = % {&GenUtils::readLocsFormat(\*ANNOT)};

my %scores = ();
my $line = <SCORE>;
while(defined($line) && $line  =~ />(\S+)\s(\d+)/) {
    my $cid = $1;
    my $clen = $2;
    my $strand = STRAND_FWD;
    
    while(defined($line = <SCORE>) && $line  !~ />(\S+)\s(\d+)/) {
        my $sigPos;
        my $sigScore;
        if($line =~ /^(\d+)\s+([-,\d,.]+)\s*$/) {
            $sigPos = $1;
            $sigScore = $2;
            if($strand == STRAND_FWD) {
		$scores{$cid}->{$sigPos} = $sigScore;
#		print "$cid $sigPos $sigScore\n";
            }
            else {
		$scores{$cid}->{$clen - $sigPos + 1} = $sigScore;
#		print "$cid $sigPos $sigScore\n";
            }
        }
        elsif($line =~ /\/\//) {
            $strand = STRAND_COMP;
        }
    }
}

close(SCORE);

my %starts = ();
foreach my $id (keys %pgenes) {
    my $gene = $pgenes{$id};
    my $cid = $pgenes{$id}->{'Exons'}->[0]->[0];
    my $start = ($gene->{'Strand'} == STRAND_FWD ?
		 $gene->{'Exons'}->[0]->[1] :
		 $gene->{'Exons'}->[-1]->[1]);
    $starts{$cid}->{$start} = $gene;
}

foreach my $id (keys %genes) {
    my $gene = $genes{$id};
    my $cid = $genes{$id}->{'Exons'}->[0]->[0];
    my $start = ($gene->{'Strand'} == STRAND_FWD ?
		 $gene->{'Exons'}->[0]->[1] :
		 $gene->{'Exons'}->[-1]->[1]);

    if(!defined($starts{$cid}->{$start})) {
 	if(!defined($scores{$cid}->{$start})) {
	    print "(0)\t\t";
	}
	else {
	    print "(", $scores{$cid}->{$start}, ")\t\t";
	}

	$gene->displayHeader(\*STDOUT);

	foreach my $id2 (keys %{$starts{$cid}}) {
	    print "\n";

	    if(!defined($scores{$cid}->{$id2})) {
		print "(0)\t\t";
	    }
	    else {
		print "(", $scores{$cid}->{$id2}, ")\t\t";
	    }
	    $starts{$cid}->{$id2}->displayHeader(\*STDOUT);
	}
	print "\n\n";

    }
}
