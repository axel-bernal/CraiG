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

my $usage = "usage : $0 subcontiglocs < locs\n";

my $contigFile = shift || die $usage;

open(SUBCONTIGS, "<$contigFile") || die "file $contigFile not found\n";
my %subcontigs = % { &GenUtils::readLocsFormat(\*SUBCONTIGS, 0, 'TRAIN', undef)};
my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};


foreach my $id (keys %genes) {
    $genes{$id}->moveLocation2OrigContigs(\%subcontigs);
#    for(my $i = 0; $i < scalar(@$exons); $i++) {
#	$exons->[$i]->[0] =  $contig->[0];
#	$exons->[$i]->[1] +=  $contig->[1] - 1;
#	$exons->[$i]->[2] +=  $contig->[1] - 1;
#    }
    $genes{$id}->dump(\*STDOUT, undef, 1);
}
