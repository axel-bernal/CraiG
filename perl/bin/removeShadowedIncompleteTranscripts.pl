#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;
use GenUtils;

my $usage = "usage : $0 shadowed_genes < shadowing_genes\n";
my $shwedFN = shift (@ARGV) || die $usage;
my %shwingGenes = % {&GenUtils::readLocsFormat(\*STDIN)};

open(GENES, $shwedFN) || die "Can't open $shwedFN\n";
my %shwedGenes = % {&GenUtils::readLocsFormat(\*GENES)};
close(GENES);

my $gene;
my $gene2;
my %discarded = ();
foreach my $id (keys %shwingGenes) {
    my $dumped = 0;
    $gene = $shwingGenes{$id};
    foreach my $id2 (keys %shwedGenes) {
        $gene2 = $shwedGenes{$id2};
        
        if($gene->contains($gene2)) {
            $discarded{$id2} = 1;
#            print STDERR "ok!\n";
            if(!$dumped) {
                $dumped = 1;
                $gene->dump(\*STDOUT, undef,  1);
            }
        }
        else {
#            print STDERR "no ok!\n";
        }
    }
}

foreach my $id2 (keys %shwedGenes) {
    $gene2 = $shwedGenes{$id2};
    if(!$discarded{$id2}) {
        $gene2->dump(\*STDOUT, undef,  1);
    }
}
