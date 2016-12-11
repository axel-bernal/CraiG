#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 [-l] < subcontigs > new_subcontigs\n";


my $subContigsAreLocs = 0;

&GetOptions("-l" => \$subContigsAreLocs,
            );

my $scLengths;
my $contigRanges;

if($subContigsAreLocs) {
    ($contigRanges, $scLengths) = @{&GenUtils::loadSubContigLocs(\*STDIN)};
}
else {
    ($contigRanges, $scLengths) = @{&GenUtils::loadSubContigSeqs(\*STDIN)};    
}

my $cid;
my $index;
my $num_subcs;

foreach $cid (keys %$contigRanges) {
    my $scs = $contigRanges->{$cid};
    $num_subcs =  scalar(@{$scs});

    if($num_subcs <= 1) {
	next;
    }

    my $lsubc = $scs->[0];
    my @marked4print = ();

    for($index = 0; $index < $num_subcs; $index++) {
	push(@marked4print, 0);
    }

    for($index = 1; $index < $num_subcs; $index++) {
	my $rsubc = $scs->[$index];
	my $olap = $lsubc->[3] - $rsubc->[2] + 1;

	if($olap > 0) {
	    $contigRanges->{$cid}->[$index - 1]->[3] = $rsubc->[2] - 1;
	    $marked4print[$index - 1] = 1;
	    $marked4print[$index] = 1;
 	}
	$lsubc = $rsubc;
    }

    for($index = 0; $index < $num_subcs; $index++) {
	if($marked4print[$index]) {
	    my $range = $contigRanges->{$cid}->[$index];
	    print ">$range->[0]\t$cid", "_", 
	    $range->[2], "_",
	    $range->[3],"\n";
	}
    }
}


#foreach $cid (keys %nsubcs) {
#    $num_subcs = scalar(@{$nsubcs{$cid}});
#    for($index = 0; $index < $num_subcs; $index++) {
#	if($nsubcs{$cid}->[$index] >= 0) {
#	    print ">$nsubcs{$cid}->[$index]->[0]\t$cid","_", $nsubcs{$cid}->[$index]->[2],"_",$nsubcs{$cid}->[$index]->[3],"\n";
#	}
#    }
#}
