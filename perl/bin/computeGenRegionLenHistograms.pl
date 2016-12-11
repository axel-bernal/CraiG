#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";

use lib "$ENV{CRAIG_HOME}/perl/lib";
use GenUtils;
# -*- perl -*-

$usage = "usage: $0 -f hist_file < training_set_of_sequences. creates hist_file.exlen his_file.inlen and hist_file.exnum smoothed histograms for the length and the number of segments appearing or each fasta entry\n";

if(defined($ARGV[0]) && ($ARGV[0] eq "-f")) {
    shift; $hist_file = shift || die $usage; 
    open(EXLEN_HIST, ">$hist_file.exlen") || die  "Cannot open $hist_file.exlen\n";
    open(INLEN_HIST, ">$hist_file.inlen") || die  "Cannot open $hist_file.inlen\n";
    open(NUM_HIST, ">$hist_file.exnum") || die  "Cannot open $hist_file.exnum\n";
}
else {
    die $usage;
}

my $lenTotal = 0;
my $numTotal = 0;
my %exlenVals = ();
my %numVals = ();
my %inlenVals = ();
my %genes = % {&GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};

for my $id (keys %genes) {
    $gene = $genes{$id};
    my $lastEnd;
    my $exons = $gene->{'Exons'};
    my ($finalExon, $initExon);
    
    $initExon = 0;
    $finalExon = scalar(@$exons);
    $lastEnd = -1;
    
    for($i = $initExon; $i < $finalExon;$i++) {
	$exlenVals{abs($exons->[$i]->[2] - $exons->[$i]->[1]) + 1}++;
	if($lastEnd > 0) {
	    $inlenVals{abs($exons->[$i]->[1] - $lastEnd)}++;
	}
	$lastEnd = $exons->[$i]->[2];
    }
    $numVals{scalar(@$exons)}++;
    
}

foreach $value (sort {$a <=> $b} keys %numVals) {
    print NUM_HIST $value," ", $numVals{$value},"\n";
}

foreach $value (sort {$a <=> $b} keys %exlenVals) {
    print EXLEN_HIST $value," ", $exlenVals{$value},"\n";
}

foreach $value (sort {$a <=> $b} keys %inlenVals) {
    print INLEN_HIST $value," ", $inlenVals{$value},"\n";
}


