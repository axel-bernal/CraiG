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
my $igenicsFile = "";
my $maxcLen = 1000000000;
my $outC;
my $contigs = "";
my $GCClasses = "43,51,57";
my $geoMeans = "";
my @geoMeans = (50000, 12000, 6000, 2500);
my $usage = "usage : $0 -ctg contigs  -GC gcClasses -means geoMeans -max max_contig_length -igenics igenicSubLocs -outC newContigSubLocs < geneLocations > newGeneLocations\n";

&GetOptions("-outC=s" => \$outC,
	    "-means=s" => \$geoMeans,
	    "-ctg=s" => \$contigs,
	    "-GC=s" => \$GCClasses,
            "-igenics=s" => \$igenicsFile,
            "-max=i" => \$maxcLen,
            );

($igenicsFile ne "" && $contigs ne "") || die "igenics and contigs must be specified. ".$usage;

if(length($geoMeans)) {
    @geoMeans = split(/\,/, $geoMeans);
}

my $i;

my @GCClasses = ();
my $lastGC = 0;
foreach my $gc (split(/\,/, $GCClasses)) {
    push(@GCClasses, [$lastGC, $gc]);
    $lastGC = $gc;
}
push(@GCClasses, [$lastGC, 100]);
($#GCClasses == $#geoMeans) || die;

my @accProbs = ();
for($i = 0; $i <= $#geoMeans; $i++) {
    push(@accProbs, \@{&Utils::getGeoAccumProbs($geoMeans[$i], $geoMeans[$i]*6)});
}

open(IGENICS, "$igenicsFile") || die "Can't open $igenicsFile\n";
open(CONTIGS, "$contigs") || die "Can't open $contigs\n";
open(OUTC, ">$outC") || die "Can't open $outC\n";
my $org = Organism->new();
my %contigs = % {$org->loadContigs(\*CONTIGS)};
close(CONTIGS);

my %igenics = %{&GenUtils::readLocsFormat(\*IGENICS, 0, 'IGENICS', undef)};
my %genes = %{&GenUtils::readLocsFormat(\*STDIN, 0, 'GENES', undef)};
my %outContigs = ();

my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);
my @igenics = sort {&Gene::compareTo($a, $b);} (values %igenics);

if(!scalar(@igenics)) {
    return;
}

my $j = 0;
my $igenicL = $igenics[0];
my $lastIgenicEnd = $igenicL->lowestCoord();

for($i = 1; $i < scalar(@igenics); $i++) {
    my $cidL =  $igenicL->contigId();
    my $igenicR = $igenics[$i];
 
    if($igenicR->contigId() ne $cidL) {
        $igenicL = $igenicR;
	$lastIgenicEnd = $igenicL->lowestCoord();
        next;
    }

    my $c = $igenicL->contigId();
#    print STDERR "$cidL\n";
    my $seq = $contigs{$c}->getSequence(0, $igenicL->highestCoord(), $igenicR->lowestCoord(), 0, 0, STRAND_FWD);
    my $gcPerc = SigUtils::getGCPerc($seq);
    my $GC = 0;
    for( ; $GC <= $#GCClasses; $GC++) {
	if($gcPerc >= $GCClasses[$GC]->[0] && $gcPerc < $GCClasses[$GC]->[1]) {
	    last;
	}
    }
    
    # finding the intergenic lengths, starting with the smaller mean dist.
    my $cbeg = 0;
    my $cend = 0;

    my $len3 = Utils::pickAtRandom($accProbs[$GC]) + 30;
    $cbeg = $igenicL->highestCoord() - $len3;
    
    if($cbeg < $lastIgenicEnd) {
	$cbeg = ($lastIgenicEnd ? $lastIgenicEnd : 1);
    }
    
    my $len5 = Utils::pickAtRandom($accProbs[$GC]) + 30;
    $cend = $igenicR->lowestCoord() + $len5;
    
    if($cend > $igenicR->highestCoord() - 1) {
	$cend = $igenicR->highestCoord() > $contigs{$c}->{'Len'} ?
	    $contigs{$c}->{'Len'} : $igenicR->highestCoord();
    }
    
#    my $cbeg = $igenicL->lowestCoord();
#    my $cend = $igenicR->highestCoord();
    
    if($cend - $cbeg + 1 > $maxcLen) {
        # try moving cbeg and end
        my $mincLen = $igenicR->lowestCoord() - $igenicL->highestCoord();
#        print STDERR "\nlarge $igenicL->{'Id'} [", $igenicL->lowestCoord()," ", $igenicL->highestCoord(), "] [", $igenicR->lowestCoord()," ", $igenicR->highestCoord(), "] [", "$cbeg $cend $mincLen]\n";

        my $toComplete = $maxcLen - $mincLen;
        if($toComplete > 0) {
            my $lenIgenicL = $igenicL->highestCoord() - $igenicL->lowestCoord();
            my $lenIgenicR = $igenicR->highestCoord() - $igenicR->lowestCoord();
            $cbeg = $igenicL->highestCoord() -
                int((1.0*$lenIgenicL/($lenIgenicL + $lenIgenicR))*$toComplete);
            $cend = $igenicR->lowestCoord() +
                int((1.0*$lenIgenicR/($lenIgenicL + $lenIgenicR))*$toComplete);
#            print STDERR "yes [$cbeg $cend]\n";
        }
        else {
            $igenicL = $igenicR;
	    $lastIgenicEnd = $cend + 1;
            next;
        }
    }

    my $new_cid = "";
#    print STDERR "igenic $igenicL->{'Id'} $cbeg $cend\n";

    while($j < scalar(@genes)) {
        my $gene = $genes[$j];
        my $cid = $gene->contigId();
        
        defined($cid) || die "undefined id for $gene->{'Id'}\n";

        my $comp = $cid cmp $cidL;
        if($comp < 0) {
            $j++;
            next;
        }
        elsif($comp > 0) {
            last;
        }

        my $glowest = $gene->lowestCoord();
        my $ghighest = $gene->highestCoord();
#        print STDERR "\tgene $gene->{'Id'} $glowest $ghighest\n";
        if($glowest < $cbeg) {
            warn "Error: Gene $gene->{'Id'} [$glowest $ghighest] overlaps intergenic regions [$cbeg $cend](removed)\n";
            $j++;
            next;
        }
        
        if($ghighest <= $cend) {
            $new_cid = $cid;
            my $offset = 1;
            if($cid =~ /^([^\.]+)\.(\d+)$/) {
                $new_cid = $1;
                $offset = $2;
            }
            $new_cid = $new_cid.".".($offset + $cbeg - 1);

            $gene->shiftCoords(-$cbeg + 1, $cend - $cbeg + 1, $new_cid);
            $gene->dump(\*STDOUT, undef, 1);
        }
        else {
            last;
        }
        $j++;
    }
    
    if(length($new_cid)) {
        print OUTC ">$new_cid\t$cidL", "_", $cbeg, "_", "$cend\n";
    }

    $igenicL = $igenicR;
    $lastIgenicEnd = $cend + 1;
}


close(IGENICS);
close(OUTC);
