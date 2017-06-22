#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $usage = "usage : $0 < contigFile\n";


my $org = Organism->new();

my %contigs = % {$org->loadContigs(\*STDIN)};
my $maskTotal = 0;
my $total = 0;
foreach my $c (keys %contigs) {
    my $mask = 0;
    my  $length = length($contigs{$c}->{'Seq'});
    for(my $i = 0; $i < $length; $i++) {
	my $cp = substr($contigs{$c}->{'Seq'}, $i, 1);
	if($cp =~ /[a-z]/) {
	    $mask++;
	}
    }
    print $c, "\t", 1.0*$mask/$length, "\n";
    $maskTotal += $mask;
    $total += $length;
}
print "Total = ", 1.0*$maskTotal/$total, "\n";
