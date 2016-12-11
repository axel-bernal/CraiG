#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $usage = "usage : $0 contigs < subcontig_info > sub_contigs\n";

my $contigFile = shift|| die $usage;
my $org = Organism->new();
my $org2 = Organism->new();
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my %subcontigs = % {$org2->loadContigs(\*STDIN)};
foreach my $id (keys %subcontigs) {
    $id =~ /^(\S+)\.(\d+)$/;
    my $len = $subcontigs{$id}->{'Len'};
    my $seq = $contigs{$1}->getSequence(0, $2, $2 + $len - 1, 0, 1, STRAND_FWD);
    print ">$id\n", $seq, "\n";
}

