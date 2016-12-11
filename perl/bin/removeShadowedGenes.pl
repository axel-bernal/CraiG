#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $gene;
my $lgene;
my $xalt = 0;
my $usage = "usage : $0  [-xalt] < locations\n-xalt excludes alternative transcripts from being removed\n";

&GetOptions("-xalt" => \$xalt);

my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);



if(!scalar(@genes)) {
    exit;
}

my $i = 0;
$lgene = $genes[$i++];
while($i < scalar(@genes)) {
    $gene = $genes[$i];

    my $b = $lgene->compareTo($gene);

    if($gene->contigId() eq $lgene->contigId()) {
        if($lgene->highestCoord() < $gene->lowestCoord() ||
           $gene->highestCoord() < $lgene->lowestCoord()) {
            ;
        }
        else {
            my $remove = 0;
            if($xalt) {
                if($gene->{'gId'} eq $lgene->{'gId'}) {
                    ;
                }
                else {
                    $remove = 1;
                }
            }
            else {
                $remove = 1;
            }

            if($remove) {
                print "overlap $lgene->{'Id'} $gene->{'Id'}\n";
                if($gene->highestCoord() - $gene->lowestCoord() >
                   $lgene->highestCoord() - $lgene->lowestCoord()) {
                    splice(@genes, $i - 1, 1);
                    $lgene = $gene;
                }
                else {
                    splice(@genes, $i, 1);
                }
                next;
            }
        }
    }
    $lgene = $gene;
    $i++;
}


foreach $gene (@genes) {
    $gene->dump(\*STDOUT, undef, 1);
}
