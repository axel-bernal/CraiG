#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 [-l] -ctg subcontig_info < score_filter > new_score_filter\nA faster program that requires much less memory to run";

my $subcontig = "";
my $subContigsAreLocs = 0;

&GetOptions("-l" => \$subContigsAreLocs,
            "-ctg=s" => \$subcontig,
            );

open(SUBCONTIG, "$subcontig") ||  die $usage;

my $scLengths;
my $contigRanges;
my %seen = ();

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
    
    if(!defined($contigRanges->{$cid})) {
        while(defined($line = <STDIN>) && $line  !~ />(\S+)\s(\d+)/) { }
        next;
    }
    
    
    my $index = 0;
    my $range = $contigRanges->{$cid}->[$index];
    my $length =  $range->[3] - $range->[2] + 1;
    my $scores;
    
    while(defined($line = <STDIN>) && $line  !~ />(\S+)\s(\d+)/) {
        if(!defined($seen{$range->[0]})) {
            print ">$range->[0]\t$length\n";
            $seen{$range->[0]} = 1;
        }

        my $sigPos;
        my $sigScore;
        
        if($line =~ /^(\d+)\s+([-,\d,.]+)\s*$/) {
            $sigPos = $1;
            $sigScore = $2;
            if($sigPos < $range->[2]) {
                
            }
            elsif($sigPos <= $range->[3]) {
                if($sigScore != 0) {
                    print $sigPos - $range->[2] + 1, "\t", $sigScore, "\n";
                }
            }
            else {
                if($index < scalar(@{$contigRanges->{$cid}}) - 1) {
                    $index++;
                    $range = $contigRanges->{$cid}->[$index];
                    $length =  $range->[3] - $range->[2] + 1;
                }
            }
        }
        elsif($line =~ /\/\//) {
            die "Only forward strand filters allowed\n";
        }
    }
}

