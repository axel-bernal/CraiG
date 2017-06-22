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
my $minSize = 0;
my $outC;
my $inC;
my $locs = 0;
my $usage = "usage : $0 -l -min minSize -outC newContig -inC contigFile < locations > newLocations\nSplits a region containing genes in multiple regions, each of them containing at most 1 gene. If -l is specified, newContig will have the new contigs as locations.\n";

&GetOptions("-outC=s" => \$outC,
            "-inC=s" => \$contigFile,
            "-min=i" => \$minSize,
            "-l" => \$locs
            );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(OUT_CONTIG, ">$outC") || die "file $outC not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};

for my $c (keys %contigs) {
    my $array = &GenUtils::markIGenics($contigs{$c}, 0, 1, 'TRAIN');
    my @genes = sort {&Gene::compareTo($a, $b);} (values %{ $contigs{$c}->{'Features'}->{'TRAIN'}});

    my $gene;
    my $i;
    my $interBeg = -1;
    my $interEnd = -1;
    my $lastSplit = 1;
    my $newSplit = 1;
    my $lastGene = 0;

    my $new_c = $c;
    my $offset = 1;
    if($c =~ /^([^\.]+)\.(\d+)$/) {
        $new_c = $1;
        $offset = $2;
    }

    for($i = 1; $i <= $contigs{$c}->{'Len'}; $i++) {
        if($array->[$i] == 0 && $interBeg == -1) {
            $interBeg = $i;
        }
        if($i > 1) {
	    if($i == $contigs{$c}->{'Len'}) {
		$interEnd = $i;
	    }
            elsif($array->[$i] == 1 && $array->[$i - 1] == 0) { 
		$interEnd = $i - 1;
	    }
	}

        if($interBeg != -1 && $interEnd != -1) {
#            print ">$interBeg $interEnd\n";
            $newSplit = ($interEnd == $contigs{$c}->{'Len'}) ? $contigs{$c}->{'Len'} : int(($interBeg + $interEnd)/2);
            if($newSplit != $contigs{$c}->{'Len'} && $newSplit - $lastSplit < $minSize) {
                $interBeg = -1;
                $interEnd = -1;
                next;
            }
            my $j = $lastGene;
            my $new_cid = $new_c.".".($offset + $lastSplit - 1);

            for(; $j < scalar(@genes); $j++) {
                $gene =  $genes[$j];
                my $min = $gene->lowestCoord();
                my $max = $gene->highestCoord();
                ($min >= $lastSplit) or die "$gene->{'Id'}\n";

                if($max <= $interBeg) {
                    $gene->shiftCoords(-$lastSplit + 1, $newSplit - $lastSplit + 1, $new_cid);
                    $gene->dump(\*STDOUT, undef, 1);
                }
                else {
                    last;
                }
            }
            if($j != $lastGene) {
                if($locs) {
                    print OUT_CONTIG ">$new_cid\t$new_c",
                    "_", ($offset + $lastSplit - 1),"_", ($offset + $newSplit - 1), "\n";
                }
                else {
                    my $cString = $contigs{$c}->getSequence(0, $lastSplit, $newSplit, 0, 0, STRAND_FWD);
                    print OUT_CONTIG ">$new_cid\n$cString\n";
                }
                
                $lastSplit = $newSplit + 1;
                $lastGene = $j;
            }
            $interBeg = -1;
            $interEnd = -1;
        }
    }
}
