#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 score_filter1 score_filter2 > new_score_filter\nFilters must be ordered with respect to contig id\n";

my $ofh = select STDERR;
$| = 1;
select $ofh;

my $score_file1 = shift || die $usage;
my $score_file2 = shift || die $usage;
open(FD1, "<$score_file1") || die $usage;
open(FD2, "<$score_file2") || die $usage;
my $filterVals;
my %filterVals1 = ();
my %filterVals2 = ();
my $line1 = <FD1>;
my $line2 = <FD2>;
my $cid1 = GenUtils::loadFilterVals4Id(\*FD1, \$line1, \%filterVals1, 0);
my $cid2 = GenUtils::loadFilterVals4Id(\*FD2, \$line2, \%filterVals2, 0);
my $strand;
my $clen;

while(1) {    
    my %filterVals = ();

    if($cid1 eq "" && $cid2 eq "") {
	last;
    }
    elsif($cid1 eq "") {
	GenUtils::printScoreFilterVals(\%filterVals2, \*STDOUT, STRAND_COMP);
	%filterVals2 = ();
	$cid2 = loadFilterVals4Id(\*FD2, \$line2, \%filterVals2, 0);
    }
    elsif($cid2 eq "") {
	GenUtils::printScoreFilterVals(\%filterVals1, \*STDOUT, STRAND_COMP);
	%filterVals1 = ();
	$cid1 = loadFilterVals4Id(\*FD1, \$line1, \%filterVals1, 0);
    }
    else {
	my $b = $cid1 cmp $cid2;
	if(!$b) {
	    $clen = scalar(@{$filterVals1{$cid1}->[0]}) - 1;
	    GenUtils::createFilterValArrays(\%filterVals, $cid1, $clen, 0);
	    
	    for($strand = 0; $strand < 2; $strand++) {
		my $filterVals1 = $filterVals1{$cid1}->[$strand];
		my $filterVals2 = $filterVals2{$cid2}->[$strand];
#		print STDERR "$cid1 $cid2 $clen $strand\n";
		for(my $i = 1; $i <= $clen; $i++) {
#		    print "$strand $i $filterVals1->[$i] $filterVals2->[$i]\n";
		    my $scoreVal = 0;		    
		    my $fVal1 = $filterVals1->[$i];
		    my $fVal2 = $filterVals2->[$i];
		    
		    if($fVal1 != 0 && $fVal2 != 0) {
			if($fVal1 != $fVal2) {	
			    die "scores are different $cid1\[$strand, $i\] = \"$scoreVal\" \"$filterVals2->[$i]\"\n";
			}
		    }
		    
		    if($fVal1 != 0) {
			$scoreVal = $filterVals1->[$i];
		    }
		    if($fVal2 != 0) {
			$scoreVal = $filterVals2->[$i];
		    }
		    
		    $filterVals{$cid1}->[$strand]->[$i] = $scoreVal;
		}
	    }
	    
	    GenUtils::printScoreFilterVals(\%filterVals, \*STDOUT, STRAND_COMP);
	    %filterVals1 = ();
	    $cid1 = loadFilterVals4Id(\*FD1, \$line1, \%filterVals1, 0);
	    %filterVals2 = ();
	    $cid2 = loadFilterVals4Id(\*FD2, \$line2, \%filterVals2, 0);
	}
	elsif($b < 0) {
	    GenUtils::printScoreFilterVals(\%filterVals1, \*STDOUT, STRAND_COMP);
	    %filterVals1 = ();
	    $cid1 = loadFilterVals4Id(\*FD1, \$line1, \%filterVals1, 0);
	}
	else {
	    GenUtils::printScoreFilterVals(\%filterVals2, \*STDOUT, STRAND_COMP);
	    %filterVals2 = ();
	    $cid2 = loadFilterVals4Id(\*FD2, \$line2, \%filterVals2, 0);
	}
    }    
}

print STDERR "\nDone\n";
exit 0;
