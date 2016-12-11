#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $usage = "usage : $0  < fasta file. Computes the GC content of each id in the fasta file\n";
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $org = Organism->new();
my %contigs = % {$org->loadContigs(\*STDIN)};

foreach my $c (keys %contigs) {
    print "$c -> ", SigUtils::getGCPerc($contigs{$c}->{'Seq'}),"\n";
}
