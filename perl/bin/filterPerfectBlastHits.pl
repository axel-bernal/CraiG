#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use strict;

my $threshold = 0.99;
my $usage = "usage : $0 -t threshold[0 .. 0.99] < blastFile > new_genes";

&GetOptions("-t=f" => \$threshold,
            );

my $line;
my $sbjt_identical = 0;
my $sbjt_beg = -1;
my $sbjt_end = -1;
my $q_id = "";
my $sbjt_strand = 'Plus';
my $q_strand = 'Plus';
my $subj = "";
my $length = 0;

while(defined($line = <STDIN>)) {
    my $update = 0;

    if($line =~ /^Query=\s*(\S+)\s*$/) {
	if($sbjt_identical) { # store sbjt
	    print "$q_id $subj $sbjt_beg $sbjt_end $q_strand\n";
	}
	$subj= "";
	$update = 1;
	$q_id = $1;
	$line = <STDIN>;
	$line =~ /\s*\((\d+)\s+letters/ || die $line;
	$length = $1;
    }
    elsif($line =~ />(\S+)/) {
	if($sbjt_identical) { # store sbjt
	    print "$q_id $subj $sbjt_beg $sbjt_end $q_strand\n";
	}
	$subj= $1;
	$update = 1;
    }
    elsif($line =~ /\sIdentities\s+\=\s+(\d+)\/(\d+)\s*/) {
	
	if($sbjt_identical) { # store sbjt
	    print "$q_id $subj $sbjt_beg $sbjt_end $q_strand\n";
	}

	$sbjt_identical = ($1 == $2 && $2 == $length) ? 1 : 0;
	$line = <STDIN>; 
	$line =~ /\sStrand\s+\=\s+(\S+)\s+\/\s+(\S+)/ || die $line;
	$q_strand = $1;
	$sbjt_strand = $2;
	$sbjt_beg = -1;
	$sbjt_end = -1;
    }
    elsif($sbjt_identical && $line =~ /Sbjct:\s(\d+)\s+\S+\s+(\d+)/) {
	if($sbjt_strand eq 'Plus') {
	    if($sbjt_beg < 0 || $1 < $sbjt_beg) { $sbjt_beg = $1; }
	    if($sbjt_end < 0 || $2 > $sbjt_end) { $sbjt_end = $2; }
	}
	else {
	    if($sbjt_beg < 0 || $1 > $sbjt_beg) { $sbjt_beg = $1; }
	    if($sbjt_end < 0 || $2 < $sbjt_end) { $sbjt_end = $2; }
	}
    }

    if($update) {
 	$sbjt_identical = 0;
	$sbjt_strand = 'Plus';
	$q_strand = 'Plus';
	$sbjt_beg = -1;
	$sbjt_end = -1;
    }
}
