#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use Getopt::Long;

my $col = 0;
my $twoVar = "";
my $numBins = 100;
my $min = 1000000000.0;
my $max = -1000000000.0;
my $userMin = undef;
my $userMax = undef;
my $usage = "$0 -col <col_no> -twoVar <2ndFile> < 1stFile\nGenerates a histogram to be drawn using Histogram.xla.\n\t-col[0] specifies in which column is the variable whose distribution will be\n\t plotted, counting from zero.\n\t-numBins [100] specifies how many bins to create for the histogram.\n\t-twoVar tell the program to produce a 2-variable histogram. Usually one of\n\t them is the positive and the other the negative occurrences\n";

if(!defined(@ARGV)) {
    die $usage;
}
&GetOptions("-col=i" => \$col,
            "-twovar=s" => \$twoVar,
            "-bins=i" => \$numBins,
            "-max=i" => \$userMax,
            "-min=i" => \$userMin,
            );

open(TMP, "$twoVar");

if(scalar(@ARGV) > 0) {
    $col = shift @ARGV;
}
my @cols1 = ();
my @cols2 = ();
my $line;
my @line = ();
my $i;
while(defined($line = <STDIN>)) {
    @line = split(/\s+/, $line);
    push(@cols1, $line[$col]);

    if($min > $line[$col]) {
        $min = $line[$col];
    }
    if($max < $line[$col]) {
        $max = $line[$col];
    }
}

while(defined($line = <TMP>)) {
    @line = split(/\s+/, $line);
    push(@cols2, $line[$col]);
    if($min > $line[$col]) {
        $min = $line[$col];
    }
    if($max < $line[$col]) {
        $max = $line[$col];
    }
}

if(defined($userMin) && $min < $userMin) {
    $min = $userMin;
}

if(defined($userMax) && $max > $userMax) {
    $max = $userMax;
}

my $binSize = ($max - $min)/$numBins;
my @binnedCols1 = ();
my @binnedCols2 = ();
for($i = 0; $i <= $numBins; $i++) {
    push(@binnedCols1, 0);
    push(@binnedCols2, 0);
}
my $m1 = 0;
my $m2 = 0;
foreach my $q (@cols1) {
    if($q < $min || $q > $max) {
        next;
    }
    my $bin = int(($q - $min)/$binSize);
    $binnedCols1[$bin]++;    
    if($binnedCols1[$bin] > $m1) {
        $m1 = $binnedCols1[$bin];
    }
}
foreach my $q (@cols2) {
    if($q < $min || $q > $max) {
        next;
    }
    my $bin = int(($q - $min)/$binSize);
    $binnedCols2[$bin]++;    
    if($binnedCols2[$bin] > $m2) {
        $m2 = $binnedCols2[$bin];
    }
}

print $min, "\t0";
if($#cols2 >= 0) {  print "\t0\n"; }
else { print "\n";}

for($i = 0; $i <= $#binnedCols1; $i++) {
    print $min + ($i + 1)*$binSize, "\t", $binnedCols1[$i];
    if($#cols2 >= 0) {  print "\t", $binnedCols2[$i]*$m1/$m2, "\n";}
    else { print "\n";}
    print $min + ($i + 1)*$binSize, "\t", $binnedCols1[$i];
    if($#cols2 >= 0) {  print "\t", $binnedCols2[$i]*$m1/$m2, "\n";}
    else { print "\n";}
}
