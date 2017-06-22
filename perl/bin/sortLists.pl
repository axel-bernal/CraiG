#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  

my @a = ();
my $loc;

my $usage = "sortList [-ascend] < list\nDefault is descend";
my $order = -1;
$order = 1 if (defined($ARGV[1] && $ARGV[1] eq "ascend"));

while($line = <STDIN>) {
    chomp($line);
    $line =~ /[^\_]+\s+\S+_(\d+)\_(\d+)\]?/;
    push(@a, [$line, $1, $2]);
}
my @sorted = sort {$order*(abs($a->[2] - $a->[1]) <=> abs($b->[2] - $b->[1])); } @a;
foreach $loc (@sorted) {
    print abs($loc->[2] - $loc->[1]) + 1, " ", $loc->[0], "\n";
}
