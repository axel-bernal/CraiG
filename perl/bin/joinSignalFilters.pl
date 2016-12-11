#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Contig;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 signal_filter1 signal_filter2 > new_signal_filter\nSignal files must be ordered with biosort.py first\n";

my $signal_file1 = shift || die $usage;
my $signal_file2 = shift || die $usage;
open(FD1, "<$signal_file1") || die $usage;
open(FD2, "<$signal_file2") || die $usage;

my $ofh = select STDERR;
$| = 1;
select $ofh;

my %filterVals1 = ();
my %filterVals2 = ();
my $filterVals;
my $line1 = <FD1>;
my $line2 = <FD2>;
my $cid1 = GenUtils::loadFilterVals4Id(\*FD1, \$line1, \%filterVals1, '-');
my $cid2 = GenUtils::loadFilterVals4Id(\*FD2, \$line2, \%filterVals2, '-');
my $strand;
my $clen;

sub maxSignalWeight {
    my ($fval) = @_;
    my @fvals = split(/\s+/, $fval);
    my $max = 0;
    for(my $k = 1; $k <= $#fvals; $k += 3) {
	$max = $fvals[$k+2] unless $max >= $fvals[$k+2];
    }
    return $max;
}

while(1) {    
    my %filterVals = ();

    if($cid1 eq "" && $cid2 eq "") {
	last;
    }
    elsif($cid1 eq "") {
	$clen = scalar(@{$filterVals2{$cid2}->[0]}) - 1;
	GenUtils::createFilterValArrays(\%filterVals, $cid2, $clen, '-');
	GenUtils::printSignalFilterVals(\%filterVals, \*STDOUT, STRAND_COMP);
	%filterVals2 = ();
	$cid2 = loadFilterVals4Id(\*FD2, \$line2, \%filterVals2, '-');
    }
    elsif($cid2 eq "") {
	GenUtils::printSignalFilterVals(\%filterVals1, \*STDOUT, STRAND_COMP);
	%filterVals1 = ();
	$cid1 = loadFilterVals4Id(\*FD1, \$line1, \%filterVals1, '-');
    }
    else {
	my $b = $cid1 cmp $cid2;
	if(!$b) {
	    $clen = scalar(@{$filterVals1{$cid1}->[0]}) - 1;
	    GenUtils::createFilterValArrays(\%filterVals, $cid1, $clen, '-');
	    
	    for($strand = 0; $strand < 2; $strand++) {
		my $filterVals1 = $filterVals1{$cid1}->[$strand];
		my $filterVals2 = $filterVals2{$cid2}->[$strand];
#		print STDERR "$cid1 $cid2 $clen $strand\n";
		for(my $i = 1; $i <= $clen; $i++) {
#		    print "$strand $i $filterVals1->[$i] $filterVals2->[$i]\n";
		    my $sigVal = "-";
		    $filterVals1->[$i] =~ /(\S)\s*/;
		    my $fVal1 = $1;
		    $filterVals2->[$i] =~ /(\S)\s*/;
		    my $fVal2 = $1;
		    
		    if($fVal1 ne "-" && $fVal2 ne "-") {
			if($fVal1 ne $fVal2) { 
			    warn "signals are different $cid1\[$strand, $i\] = \"$filterVals1->[$i]\" and \"$filterVals2->[$i], removing lowest\"\n";
			    my $maxw1 = maxSignalWeight($filterVals1->[$i]);
			    my $maxw2 = maxSignalWeight($filterVals2->[$i]);
			    $sigVal = ($maxw1 > $maxw2 ? $filterVals1->[$i] :
				       $filterVals2->[$i]);
			}		       
		    }
		    elsif($fVal1 ne "-") {
			$sigVal = $filterVals1->[$i];
		    }
		    elsif($fVal2 ne "-") {
			$sigVal = $filterVals2->[$i];
		    }
		    
		    $filterVals{$cid1}->[$strand]->[$i] = $sigVal;
		}
	    }

	    GenUtils::printSignalFilterVals(\%filterVals, \*STDOUT, STRAND_COMP);
	    %filterVals1 = ();
	    $cid1 = loadFilterVals4Id(\*FD1, \$line1, \%filterVals1, '-');
	    %filterVals2 = ();
	    $cid2 = loadFilterVals4Id(\*FD2, \$line2, \%filterVals2, '-');
	}
	elsif($b < 0) {
	    GenUtils::printSignalFilterVals(\%filterVals1, \*STDOUT, STRAND_COMP);
	    %filterVals1 = ();
	    $cid1 = loadFilterVals4Id(\*FD1, \$line1, \%filterVals1, '-');
	}
	else {
	    my $clen = scalar(@{$filterVals1{$cid2}->[0]}) - 1;
	    GenUtils::createFilterValArrays(\%filterVals, $cid2, $clen, '-');
	    GenUtils::printSignalFilterVals(\%filterVals, \*STDOUT, STRAND_COMP);
	    %filterVals2 = ();
	    $cid2 = loadFilterVals4Id(\*FD2, \$line2, \%filterVals2, '-');
	}
    }
}

print STDERR "\nDone\n";
exit 0;
