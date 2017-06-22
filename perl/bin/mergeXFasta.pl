#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use GenUtils;
use Contig;

my $usage = "usage : $0 codr_xfasta1 codr_xfasta2 n";

my $xcodr1 = shift ||  die $usage;
my $xcodr2 = shift ||  die $usage;
open(CODR1, "<$xcodr1") || die $usage;
open(CODR2, "<$xcodr2") || die $usage;

my %codr1 = %{&GenUtils::loadExtendedFasta(\*CODR1)};
my %codr2 = %{&GenUtils::loadExtendedFasta(\*CODR2)};

foreach my $key (keys %codr1) {
    my $pos = 1;
    print ">$key\n";
    my $crArr1 = $codr1{$key};
    for(my $i = 0; $i < scalar(@$crArr1); $i++) {
        my $cr1 = $crArr1->[$i];
        if($cr1->[0] == 1) {
            my $crArr2 = $codr2{$key};
            for(my $j = 0; $j < scalar(@$crArr2); $j++) {
                my $cr2 = $crArr2->[$j];
                print "$cr2->[0] x $cr2->[1] ", $pos, "\n"; 
                $pos += $cr2->[1];
            }
        }
        else { 
            print "$cr1->[0] x $cr1->[1] ", $pos, "\n"; 
            $pos += $cr1->[1];
        }
    }
}

close(CODR1);
close(CODR2);
