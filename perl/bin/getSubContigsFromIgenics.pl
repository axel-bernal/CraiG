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
my $maxcLen = 1000000000;

my $usage = "usage : $0 [--use-subcontiglocs_as_ids] [--use-cid-offset] [--include-truncated] -max max_contig_length < igenicSubLocs >  newContigSubLocs\nIf there are no intergenic regions at the beginning or end of the contigs, then the program assumes that entries cid_0_0 and cid_length(cid)_length(cid) are present in the input\n";

my $use_subcontiglocs_as_ids = 0;
my $use_cid_offset = 0;
&GetOptions(
    "-max=i" => \$maxcLen,
    "--use-cid-offset" => \$use_cid_offset,
    "--use-subcontiglocs-as-ids" => \$use_subcontiglocs_as_ids,
    );

my %igenics = %{&GenUtils::readLocsFormat(\*STDIN, 0, 'STDIN', undef)};
my %outContigs = ();
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

    my $new_cid = $cidL;
    my $offset = $cbeg;
    if($use_cid_offset && $cidL =~ /^([^\.]+)\.(\d+)$/) {
	$new_cid = $1;
	$offset = $2;
    }
    
    $new_cid = ($use_subcontiglocs_as_ids ? $cidL."_".$cbeg."_".$cend : $new_cid.".".($offset + $cbeg - 1));

    print STDOUT ">$new_cid\t$cidL", "_", $cbeg, "_", "$cend\n";    
    $igenicL = $igenicR;    
}
