#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use strict;

$| = 1;
my $usage = "usage: $0 list_ranges < gff3 > filtered_ggf3\n";
my @ranges = ();
my $file = shift or die $usage;

open(RANGES, $file) or die $usage;

while(<RANGES>) {
    chomp();
    $_ =~ /([^\s,\:]+)\:(\d+)\.\.(\d+)/;
    push(@ranges, [$1, $2, $3]);
}

close(RANGES);

my $inRange;
my@line;
while(defined($_ = <STDIN>)) {
    
    @line = split(/\s+/, $_);
    $inRange = 0;
    foreach my $r (@ranges) {
        if($line[0] eq $r->[0] && $line[3] >= $r->[1] && $line[4] <= $r->[2]) {
            $inRange = 1;
            last;
        }
    }

    if($inRange) {
        print $_;
    }
}
