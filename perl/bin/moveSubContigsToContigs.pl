#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use GenUtils;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

# this works better
#perl -ne '$_ =~ /^(\S+)\s+(\S+)\.(\d+)_(\d+)_(\d+)$/; print "$1\t$2_", $3+$4 -1,"_",$3+$5 -1,"\n";'

my $usage = "usage : $0 contigs < subcontigs\n";

my $contigFile = shift || die $usage;

open(CONTIGS, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % { &GenUtils::readLocsFormat(\*CONTIGS, 0, 'TRAIN', undef)};
my %subcontigs = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my $subcontig;
my $contig;
my $id;

foreach $id (keys %subcontigs) {
    $subcontig = $subcontigs{$id}->{'Exons'}->[0];
    (defined($contigs{$subcontig->[0]})) || die $subcontig;
    $contig = $contigs{$subcontig->[0]}->{'Exons'}->[0];
    $subcontig->[0] =  $contig->[0];
    $subcontig->[1] +=  $contig->[1] - 1;
    $subcontig->[2] +=  $contig->[1] - 1;
    $subcontigs{$id}->dump(\*STDOUT, undef, 1);
}
