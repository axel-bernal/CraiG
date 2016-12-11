#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use GenUtils;
use strict;
use Getopt::Long;
use File::Basename;

$| = 1;

# -*- perl -*-
my $usage = "$0 [-notest] numParts locs contigs tags [file_in_fasta ..]\n";
my $notest = 0;

($#ARGV >= 0) || die $usage;

if($ARGV[0] =~ /\-notest/) {
    $notest = 1;
    shift;
}

my $numParts = shift || $usage;
my $locs = shift || $usage;

(scalar(@ARGV) > 0) || die "Nothing to do. ".$usage;

my $org = Organism->new();

my $prefix = basename($ARGV[0]);
$prefix =~ /(\S+)\.[^\.]+/;
$prefix = dirname($ARGV[0])."/".$1;

my $locsPrefix = basename($locs);
$locsPrefix =~ /(\S+)\.[^\.]+/;
$locsPrefix = $1;

`grep ">" $ARGV[0] > $prefix.ids`;
`scramble4XValid.pl -n $numParts -i $prefix.ids -fmt fasta > $prefix.sets`;

for(my $i = 0; $i < 1; $i++) {
    if(-d "$prefix.run.$i") {
        `rm $prefix.run.$i/*`;
    }
    else {
        `mkdir $prefix.run.$i`;
    }
    for(my $j = 0; $j < $numParts; $j++) {
        my $part = ($j + $i) % $numParts;
        my $name = "train";
        if($j == 0) {
            $name = "valid";
        }
        elsif($j <= 2 && !$notest) {
            $name = "test";
        }
        `choose4XValid.pl fasta $prefix.ids $prefix.sets $part | biointers.py $locs - -by ci >> "$prefix.run.$i/$locsPrefix.$name.locs"`;

        for(my $k = 0; $k < scalar(@ARGV); $k++) {
#            if($k > 2 && $name ne "test") {
#                $name = "4train";
#            }
            my $fasta = basename($ARGV[$k]);
            $fasta =~ /(\S+)\.([^\.]+)/;
#            if($2 eq "tag") {
                `choose4XValid.pl fasta $ARGV[$k] $prefix.sets $part >> "$prefix.run.$i/$1.$name.$2"`;
#            }
        }        
    }
}
