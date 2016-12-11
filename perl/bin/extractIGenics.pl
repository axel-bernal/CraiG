#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use Gene;
use strict;

use enum qw(STRAND_FWD=0 STRAND_COMP=1);
my $contigFile;
my $minSize = 1;
my $usage = "usage : $0 -min minSize -ctg contigFile < gene_locations > igenics\n";

&GetOptions(
            "-ctg=s" => \$contigFile,
            "-min=i" => \$minSize,
            );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};

for my $c (keys %contigs) {
    my $array = &GenUtils::markIGenics($contigs{$c}, 0, 1, 'TRAIN', 0);

    my $gene;
    my $i;
    my $interBeg = -1;
    my $interEnd = -1;
    my $lastSplit = 1;
    my $newSplit = 1;
    my $lastGene = 0;
    if($array->[1] != 0) {
        print ">$c.0\t$c"."_0_0\n";
    }

    for($i = 1; $i <= $contigs{$c}->{'Len'}; $i++) {
        $interBeg = $i if !$array->[$i] && $interBeg == -1;
        if($i > 1) {
	    $interEnd = $i - 1 if $array->[$i] == 1 && !$array->[$i - 1];
	    $interEnd = $i if $i == $contigs{$c}->{'Len'} && !$array->[$i];
        }

        if($interBeg != -1 && $interEnd != -1) {
#            print ">$interBeg $interEnd\n";
            
            if($interEnd != $contigs{$c}->{'Len'} && $interEnd - $interBeg + 1 < $minSize) {
                $interBeg = -1;
                $interEnd = -1;
                next;
            }
            
#            my $cString = $contigs{$c}->getSequence(0, $interBeg, $interEnd, 0, 0, STRAND_FWD);
            print ">$c.$interBeg\t$c"."_",$interBeg,"_", $interEnd, "\n";
            $interBeg = -1;
            $interEnd = -1;
        }
    }

    if($array->[$contigs{$c}->{'Len'}] != 0) {
        print ">$c.", $contigs{$c}->{'Len'} + 1, 
        "\t$c"."_", $contigs{$c}->{'Len'} + 1, "_", $contigs{$c}->{'Len'} + 1, "\n";
    }

}
