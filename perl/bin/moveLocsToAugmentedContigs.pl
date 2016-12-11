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

my $usage = "usage : $0 xfasta_aug_contigs < locs\n";

my $contigFile = shift || die $usage;

open(CONTIGS, "<$contigFile") || die "file $contigFile not found\n";
my %marks = %{&GenUtils::loadExtendedFasta(\*CONTIGS)};
my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my $gene;
my $id;

foreach $id (keys %genes) {
    $gene = $genes{$id};
    my $mark = $marks{$gene->contigId()};
    if($mark->[0]->[0] == 0) {
        my $len = 0;
        for(my $i = 0; $i < scalar(@{$mark}); $i++) {
            $len += $mark->[$i]->[1];
        }
#        print "**$id ", $gene->contigId()", $mark->[0]->[1], $len,\n";
        $gene->shiftCoords($mark->[0]->[1], $len);
    }
    $gene->dump(\*STDOUT, undef, 1);
}
