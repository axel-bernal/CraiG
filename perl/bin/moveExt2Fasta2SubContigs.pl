#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 [-l] -ctg subcontig_info < ext2fasta > new_ext2fasta\n";

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
my $fillChar;

while(defined($line) && $line  =~ />(\S+)\s(\d+)\s(\S)/) {
    my $cid = $1;
    my $clen = $2;
    $fillChar = $3;

    if(!defined($contigRanges->{$cid})) {
        while(defined($line = <STDIN>) && $line  !~ />\S+\s\d+\s\S/) { }
        next;
    }

    my $index = 0;
    my $range = $contigRanges->{$cid}->[$index];
    my $length =  $range->[3] - $range->[2] + 1;

    if(!defined($seen{$range->[0]})) {
	print ">$range->[0]\t$length\t$fillChar\n";
	$seen{$range->[0]} = 1;
    }
    
    while(defined($line = <STDIN>) && $line  !~ />(\S+)\s(\d+)\s(\S)/) {

        my $pos;
        my $seq;
        
        if($line =~ /^(\d+)\s+(\S+)$/) {
            $pos = $1;
            $seq = $2;
	    while($pos > $range->[3] - 1 &&
		  ($index < scalar(@{$contigRanges->{$cid}}) - 1)) {
		$index++;
		$range = $contigRanges->{$cid}->[$index];
		$length =  $range->[3] - $range->[2] + 1;
		
		if(!defined($seen{$range->[0]})) {
		    print ">$range->[0]\t$length\t$fillChar\n";
		    $seen{$range->[0]} = 1;
		}		
            }

            if($pos + length($seq) <= $range->[2] - 1) {
                ;
            }
            elsif($pos < $range->[2] - 1) {
                my $thisLen = Utils::min(length($seq) - ($range->[2] - 1 - $pos),
                                         $length);
                print "0\t", substr($seq, $range->[2] - 1 - $pos, $thisLen), "\n";
            }
            elsif($pos <= $range->[3] - 1) {
                print $pos - ($range->[2] - 1), "\t", 
                substr($seq, 0, $range->[3] - $pos), "\n";
            }
        }
        elsif($line =~ /\/\//) {
            die "Only forward strand filters allowed\n";
        }
    }
}

foreach my $cid (keys %$contigRanges) {
    for(my $index = 0; $index < scalar(@{$contigRanges->{$cid}}); $index++) {
        my $range = $contigRanges->{$cid}->[$index];
        my $length =  $range->[3] - $range->[2] + 1;        
        if(!defined($seen{$range->[0]})) {
            print ">$range->[0]\t$length\t$fillChar\n";
            $seen{$range->[0]} = 1;
        }        
    }
}

