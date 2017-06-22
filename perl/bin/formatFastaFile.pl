#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < contigFile \n";

my $org = Organism->new();

my %contigs = % {$org->loadContigs(\*STDIN)};

foreach my $id (keys %contigs) {
    $contigs{$id}->{'Seq'} =~ s/\s+//g;
    print ">$id\n", $contigs{$id}->{'Seq'}, "\n";
}

