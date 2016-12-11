#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $usage = "usage : $0 contig  < locations\n";
my $contigFile = shift ||  die $usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};

for my $c (keys %contigs) {
    my $array = &GenUtils::markIGenics($contigs{$c}, 0, 1, 'TRAIN');
    my $interBeg = -1;
    my $interEnd = -1;
    for(my $i = 1; $i <= $contigs{$c}->{'Len'}; $i++) {
        if($array->[$i] == 0 && $interBeg == -1) {
            $interBeg = $i;
        }
        if($i > 1) {
            if(($array->[$i] == 1 && $array->[$i - 1] == 0) || $i == $contigs{$c}->{'Len'}) {
                $interEnd = $i - 1;
            }
        }

        if($interBeg != -1 && $interEnd != -1) {
            print $c."_".$interBeg."_$interEnd\n";
            $interBeg = -1;
            $interEnd = -1;
        }        
    }
}
