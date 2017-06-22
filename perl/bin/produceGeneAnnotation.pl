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
my $include_truncated = 0;
my $outC;
my $usage = "usage : $0 [--include-truncated] -max max_contig_length -igenics igenicSubLocs -outC newContigSubLocs < geneLocations > newGeneLocations\n";

&GetOptions("-outC=s" => \$outC,
            "-igenics=s" => \$igenicsFile,
	    "--include-truncated" => \$include_truncated,
            "-max=i" => \$maxcLen,
            );

($igenicsFile ne "") || die "igenic file must be specified. ".$usage;

open(IGENICS, "$igenicsFile") || die "Can't open $igenicsFile\n";
open(OUTC, ">$outC") || die "Can't open $outC\n";

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

for(my $i = 1; $i < scalar(@igenics); $i++) {
    my $cidL =  $igenicL->contigId();
    my $igenicR = $igenics[$i];
 
    if($igenicR->contigId() ne $cidL) {
        $igenicL = $igenicR;
        next;
    }

    my $cbeg = int(($igenicL->highestCoord() + $igenicL->lowestCoord())/2) + 1;
    my $cend = int(($igenicR->highestCoord() + $igenicR->lowestCoord())/2) - 1;
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
            next;
        }
    }

    my $new_cid = "";
#    print STDERR "igenic $igenicL->{'Id'} $cbeg $cend\n";

    while($j < scalar(@genes)) {
        my $gene = Gene->new('Gene' => $genes[$j]);
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
	my $remove = 0;
        if($glowest < $cbeg) {
	    warn "Warning: Gene $gene->{'Id'} [$glowest $ghighest] overlaps intergenic regions [$cbeg $cend]";
	    if($ghighest < $cbeg) {
		$remove = 1;
	    }
	    else {
		if($include_truncated) { # truncate gene
		    warn "(truncated)\n";
		    $gene->cropCoords($cbeg, $cend);
		    $glowest = $gene->highestCoord();
		    $ghighest = $gene->highestCoord();
		}
		else {
		    $remove = 1;
		}
	    }	    
	}

	if($remove) {
	    warn "(removed)\n";
	    $j++;
	    next;
	}
	
#	print STDERR "inserting $gene->{'Id'}($glowest, $ghighest) in [$cbeg $cend]\n";
        if($ghighest <= $cend && $ghighest > 0) {
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

}


close(IGENICS);
close(OUTC);
