#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use Getopt::Long;

my $str_bin_sizes = "1";
my $log_scale = 0;
my @bin_sizes = ();
my $usage = "$0 --bin-size COMMA_SEP_BINS[1] --log-scale < Input File > Histogram\nGenerates a histogram for the input file. The number of columns is determined by the number of BINS\n";

&GetOptions("--bin-size=s" => \$str_bin_sizes,
	    "--log-scale" => \$log_scale,
            );

@bin_sizes = split(/\,/, $str_bin_sizes);

my %histogram = ();
my $line;
my @line = ();
my @boundary_vals = ();
my $i;
my $key;

for($i = 0; $i <= $#bin_sizes; $i++) {
    push(@boundary_vals, [1000000, -1000000]);
}

while(defined($line = <STDIN>)) {
    @line = split(/\s+/, $line);
    die "number of bin sizes should match the number of columns\n" if $#line != $#bin_sizes;
    $key = "";
    for($i = 0; $i <= $#line; $i++) {
	my $oval = $log_scale ?
	    ($line[$i] ? log($line[$i])/log(10.0) : -2) :
	    $line[$i];
	my $bin = int($oval/$bin_sizes[$i]);
	my $val = $bin_sizes[$i]*$bin;
	$boundary_vals[$i]->[0] = $val if $boundary_vals[$i]->[0] > $val;
	$boundary_vals[$i]->[1] = $val if $boundary_vals[$i]->[1] < $val;
	$key .= $val;
	$key .= "\t" if $i < $#line;
    }

    $histogram{$key}++;
}

foreach $key (keys %histogram) {
    print $key,"\t", $log_scale ? log($histogram{$key}) : $histogram{$key}, "\n";
}

print STDERR "boundary vals for each column are:\n";
for($i = 0; $i <= $#bin_sizes; $i++) {
    print STDERR "column $i: ", join("\t", @{$boundary_vals[$i]}), "\n";
}
