#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use GenUtils;


my $maxSize = 800000;
my $maxIntronSize = 500000;
my $usage = "usage : $0 minGeneLen < genes\n";
my $minGeneSize = shift || die $usage;
my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
my $id;
OUT: foreach $id (keys %genes) {
    my $exons = $genes{$id}->{'Exons'};
    my $size = 0;
    if(abs($exons->[-1]->[2] - $exons->[0]->[1]) >= $maxSize) {
        next;
    }
    for(my $i = 0; $i <= $#$exons; $i++) {
        $size += abs($exons->[$i]->[2] - $exons->[$i]->[1]) + 1;
#        if($i == 0) {
#            next;
#        }
#        if(abs($exons->[$i]->[1] - $exons->[$i - 1]->[2]) > $maxIntronSize) {
#            next OUT;
#        }
    }
    if($size < $minGeneSize) {
        next;
    }
    $genes{$id}->dump(\*STDOUT, undef, 1);
}
