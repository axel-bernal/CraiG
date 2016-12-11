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

my $usage = "usage : $0  lens  < score_filter > new_score_filter\n";

my $lensFile = shift || die $usage;
my %lens = ();
open(LENGTHS, "$lensFile") ||  die $usage;
my $line;
while(defined($line = <LENGTHS>)) {
    $line =~ /(\d+)\s+(\S+)/;
    $lens{$2} = $1;
}
close(LENGTHS);

my %fVals = ();
$line = <STDIN>;
my $cid = GenUtils::loadFilterVals4Id(\*STDIN, \$line, \%fVals, 0);
my $strand;

while(1) {
    if($cid eq "") {
	last;
    }    
    my $clen = $lens{$cid};
    my $f_vals = $fVals{$cid};
    
    for(my $i = 1; $i <= $clen; $i++) {
#	    print "$strand $i $f_vals1->[$i] $f_vals2->[$i]\n";
	my $score_val = $f_vals->[STRAND_COMP]->[$i];
	$f_vals->[STRAND_COMP]->[$i] = 0;
	$f_vals->[STRAND_FWD]->[$clen - $i + 1] += $score_val;
    }	

    GenUtils::printContigScoreFilterVals($cid, \%fVals, \*STDOUT, STRAND_FWD);
    %fVals = ();
    $cid = loadFilterVals4Id(\*STDIN, \$line, \%fVals, 0);

}


