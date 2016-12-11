#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);
my $pflocs = "";
my $gene;
my $pfgene;
my $usage = "usage : $0 pfam_locations < locations > ext_locs\n\t\"extends\" locations following the five prime direction,\n\tusing pfam_locations' start as a guide\n\tIt outputs the extended locations.\n";

$pflocs = shift || die $usage;
open(PFLOCS, "$pflocs") || die $usage;

my %pfgenes = % { &GenUtils::readLocsFormat(\*PFLOCS, 0, 'OLAP', undef)};
my @pfgenes = sort {&Gene::compareTo($a, $b);} (values %pfgenes);

my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);

close(PFLOCS);

my $j = 0;
my $i = 0;
my $length;
my %pfpreds = ();

for( ; $i < scalar(@genes); $i++) {
    $gene = $genes[$i];
    my $golap = 0;
#    print "[",$gene->highestCoord()," ",  $gene->lowestCoord(),"]\n";
    # find the overlapping pfam gene
    while($j < scalar(@pfgenes)) {
        $pfgene = $pfgenes[$j];

        if(!defined($pfpreds{$pfgene->contigId()})) {
            $pfpreds{$pfgene->contigId()} = [0, 0];
        }
        $pfpreds{$pfgene->contigId()}->[$pfgene->{'Strand'}] = 1;

        my $b = $pfgene->compareTo($gene);
#        print "\t[", $pfgene->highestCoord()," ",$pfgene->lowestCoord(),"] $b\n";
        if($gene->contigId() eq $pfgene->contigId() &&
           $gene->{'Strand'} == $pfgene->{'Strand'} &&
           !$pfgene->{'NoStart'}) {
            
            if($gene->highestCoord() >= $pfgene->lowestCoord() &&
               $pfgene->highestCoord() >= $gene->lowestCoord()) {
                $golap = 1;
                $j++;
                last;
            }
        }
        if($b <= 0) {  
            $j++; 
        }
        else {  
            last; 
        }
    }
    
    if(!$golap) {  next; }
          
    my $k = 0; my $l;
    my $eolap = 0;
#    print "golap $gene->{'Id'} (", scalar(@{$gene->{'Exons'}}), ") ", $pfgene->{'Id'}, " (", scalar(@{$pfgene->{'Exons'}}), ")\n";

    my $e; my $pfe;
    # find the overlapping pfam exon
    while($k < scalar(@{$gene->{'Exons'}})) {
        $e = $gene->{'Exons'}->[$k];
        $eolap = 0;
        $l = 0;

        while($l < scalar(@{$pfgene->{'Exons'}})) {
            $pfe = $pfgene->{'Exons'}->[$l];
            if($gene->{'Strand'} == STRAND_COMP) {
                if($pfe->[1] >= $e->[2] && $pfe->[2] <= $e->[1]) {
                    $eolap = 1;
                    last;
                }
            }
            else {
                if($pfe->[1] <= $e->[2] && $pfe->[2] >= $e->[1]) {
                    $eolap = 1;
                    last;
                }
            }
            $l++;
        }

        if($eolap) {  last; }
        $k++;

    }
    
    if(!$eolap) {  next; }

    # $k and $l contain overlapping exons
#    print "eolap ", join("\t", @$e), "($k) and ", join("\t", @$pfe), "($l)\n";

    my @exons = (); my $m = 0;
#    print "added ";
    for( ; $m <= $l; $m++)  {
        push(@exons, $pfgene->{'Exons'}->[$m]);
#        print join("\t", @{$pfgene->{'Exons'}->[$m]}), ", ";
    }
#    print "\nchanged $exons[-1]->[2] to ", $gene->{'Exons'}->[$k]->[2], "($k)\n";
    $exons[-1]->[2] = $gene->{'Exons'}->[$k]->[2];

#    print "remaining ";
    for($m = $k + 1; $m < scalar(@{$gene->{'Exons'}}); $m++) {
        push(@exons, $gene->{'Exons'}->[$m]);
#        print join("\t", @{$gene->{'Exons'}->[$m]}), ", ";
    } 

    $gene->{'Exons'} = \@exons;
}

foreach $gene (@genes) {
    my $notafp = 1;
    if(defined($pfpreds{$gene->contigId()}) &&
       !$pfpreds{$gene->contigId()}->[$gene->{'Strand'}] &&
       $pfpreds{$gene->contigId()}->[!$gene->{'Strand'}]) {
#        $notafp = 0;
    }
    if($notafp) {
        $gene->dump(\*STDOUT, undef, 1);
    }
}
