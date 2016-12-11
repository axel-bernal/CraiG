#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use GenUtils;
use Getopt::Long;
use strict;

my $id = "";
my $beg = 0;
my $end = 0;
my $usage = "usage : scoreRegion.pl -sc scoreFile < fasta\n";
my %scores = ();
my $scFile = "";

&GetOptions("-sc=s" => \$scFile,
	    );
open(SC,"<$scFile") || die "can't open $scFile\n$usage";

my $line = <SC>;
while(defined($line) && $line =~ /^>(\S+)/) {
    $id = $1;
    push(@{$scores{$id}}, 0);
    while(defined($line = <SC>) && $line !~ /^>/) {
	chomp($line);
	push(@{$scores{$id}}, $line); 
    }
}
close(SC);

while(defined($line = <>) && $line =~ /^>(\S+)\s+(.+)/) {
    $id = $1;
    my $i;
    my $j;
    my $k;
    my (undef, undef, $exons, undef) = @{&GenUtils::readContigLocation($2)};
    my $lastEnd = 0;
#    print "\"$2\" ", scalar(@$exons),"\n";
    if(!defined($scores{$id})) {
	next;
    }

    my @gscore1 = (0,0,0);
    my @gscore2 = (0,0,0);

    for($j = 0; $j < scalar(@$exons); $j++) {
#	print join(" ",@{$exons->[$j]}),"\n";
	my @score = (0,0,0);

	$beg = $lastEnd + 1;
	$lastEnd = $end = $beg + $exons->[$j]->[2] - $exons->[$j]->[1];
	my $frame = ($beg - 1) % 3;
	$end -= 3 if $j == scalar(@$exons) - 1;
	for($i = $beg; $i < $end; $i += 3) {
	    for($k = 0; $k < 3; $k++) {
		$score[($k + $frame)%3] += $scores{$id}->[$i+$k];
	    }
	    for($k = 0; $k < 3; $k++) {
		print $i+$k,"\t",$scores{$id}->[$i+$k],"\n";
	    }
	}
	for($k = 0; $k < 3; $k++) {
	    $gscore1[$k] += $score[$k];
	}
	print "$id\t$beg\t$end\t",join("\t", @score),"\n";
    }
    for($i = 1; $i < $lastEnd - 3; $i += 3) {
	for(my $k = 0; $k < 3; $k++) {
	    $gscore2[$k] += $scores{$id}->[$i+$k];
	}
    }
    print "$id\t1\t$lastEnd\t",join("\t", @gscore1),"\n";
    print "$id\t1\t$lastEnd\t",join("\t", @gscore2),"\n";
    $lastEnd == scalar(@{$scores{$id}}) - 1 || die "$lastEnd".(scalar(@{$scores{$id}}) - 1)."\n";
}


