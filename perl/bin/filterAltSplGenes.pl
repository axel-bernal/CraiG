#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use GenUtils;

my $usage = "usage : $0 < genes\n";

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
my %altSplices = % {&GenUtils::produceAltSplices(\%genes)};


foreach my $id (keys %altSplices) {
    if(scalar(@{$altSplices{$id}}) == 1) {
        foreach my $spl (@{$altSplices{$id}}) {
#        my $spl = $altSplices{$id}->[0];
            if($spl ne "") {
                $spl = "|".$spl;
            }
            defined($genes{$id.$spl}) or die $id.$spl."\n";
            $genes{$id.$spl}->dump(\*STDOUT, undef,  1);
        }
    }
}
