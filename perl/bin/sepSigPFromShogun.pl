#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Getopt::Long;
use strict;

my $usage = "usage : $0 prefix seqIds < sigScoreFile\n";
my $counter = 1;

my $prefix = shift @ARGV || die $usage;
my $seqIdFile = shift @ARGV ||  die $usage;
open(SEQIDS, "<$seqIdFile") || die $usage;

my $seqLine = <SEQIDS>;

defined($seqLine) || die "error $seqLine\n";
open(OUT, ">$prefix.$counter");

$seqLine =~ /(\S+)\s+(\d+)/;
my $seqId = $1;
my $seqLen = $2;

print OUT ">$seqId\t$seqLen\n";

my $offset = 50;
my $line=  <STDIN>;

while(defined($line) && $line =~ /^\S\s+(\d+)\s+([-,\d,.]+)\s*$/) {
    my $sigPos = $1 - $offset;
    my $sigScore = $2;

    if($sigPos <= $seqLen) {
        print OUT $sigPos, "\t", $sigScore, "\n";
        $line = <STDIN>;
    }
    else {
        my $seqLine = <SEQIDS>;
        
        defined($seqLine) || die "error $seqId $seqLen\n";       

        if($seqLine =~ /(\S+)\s+(\d+)/) {
            $offset = $offset + $seqLen + 50;
            $seqId = $1;
            $seqLen = $2;
            print OUT ">$seqId\t$seqLen\n";
            $line = <STDIN>;
        }
        elsif($seqLine =~ /\/\//) {
            close(OUT);
            $counter++;
            open(OUT, ">$prefix.$counter");
        }
        else {  die "error $seqId $seqLen\n"; }
        
            
    }    
}

close(OUT);
close(SEQIDS);
