#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < locations\n";

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
foreach my $id (keys %genes) {
    my $exons = $genes{$id}->{'3UTRExons'};
    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
        print "(3UTR)",$exons->[$i]->[0], "_", $exons->[$i]->[1], "_", $exons->[$i]->[2], "\n";	
    }

    $exons = $genes{$id}->{'Exons'};
    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
        if(scalar(@$exons) == 1) {  print "(single)"; }
        else {
            if($i == 0) {  print "(initial)"; }
            elsif($i == scalar(@$exons) - 1)  {  print "(terminal)"; }
            else  {  print "(internal)"; }
        }                                           
        print $exons->[$i]->[0], "_", $exons->[$i]->[1], "_", $exons->[$i]->[2], "\n";
    }

    $exons = $genes{$id}->{'5UTRExons'};
    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
        print "(5UTR)",$exons->[$i]->[0], "_", $exons->[$i]->[1], "_", $exons->[$i]->[2], "\n";	
    }

}
