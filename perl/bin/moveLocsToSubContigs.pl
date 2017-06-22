#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use GenUtils;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $usage = "usage : $0  [--preserve-ids] [-l] -ctg subcontig_info < locs\n\tUse --preserve-ids when you don't want the genes to have their\n\tsubcontig name as part of the id";
my $subContigsAreLocs = 0;
my $preserve_ids = 0;
my $contigFile = "";
&GetOptions("-l" => \$subContigsAreLocs,
            "-ctg=s" => \$contigFile,
	    "-preserve-ids" => \$preserve_ids,
            );


(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();
my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my $index;
my $gene;
my $id;
my $cid;
my $len;
my $offset;
my $exon;
my $i;

open(CONTIGS, "<$contigFile") || die "file $contigFile not found\n";

my $scLengths;
my $contigRanges;

if($subContigsAreLocs) {
    my $line;
    ($contigRanges, $scLengths) = @{&GenUtils::loadSubContigLocs(\*CONTIGS)};
}
else {
    ($contigRanges, $scLengths) = @{&GenUtils::loadSubContigSeqs(\*CONTIGS)};    
}

close(CONTIGS);

my @genes = sort {Gene::compareTo($a, $b);} values %genes;
my $last_cid = "";
#print STDERR scalar(@genes), "\n";

for($i = 0; $i < scalar(@genes); $i++) {
    $gene = $genes[$i];
    $id = $gene->{'Id'};
    $cid = $gene->contigId();
    my $sub_cid;

    if(!defined($contigRanges->{$cid})) {
        next;
    }

    if($last_cid eq "" || $last_cid ne $cid) {
        $index = 0;
    }

    my $lowestC = $gene->lowestCoord();
    my $highestC = $gene->highestCoord();
#    print "$gene->{'Id'} [", $lowestC, " ", $highestC, "]\n";
    
    while($index < scalar(@{$contigRanges->{$cid}})) {
        
        my $range = $contigRanges->{$cid}->[$index];
#        print "\t", $range->[2], " ", $range->[3], "\n";
        if($range->[2] > $highestC) {
            last;
        }

        if($range->[3] < $lowestC) {
            $index++;
            next;
        }

        my $thisGene = Gene->new('Gene' => $gene);
#        if($thisGene->lowestCoord() >= $contigRanges->{$cid}->[$index]->[1]
#           && $thisGene->highestCoord() <= $contigRanges->{$cid}->[$index]->[2])  {

	$sub_cid = $range->[0];
        my $scLength = $scLengths->{$sub_cid};
        $offset = $range->[1] - $range->[2] + 1;
        
#        print "\t", scalar(@{$thisGene->{'Exons'}}), "exons. Ok for $cid to be replaced by $sub_cid for $thisGene->{'Id'}!\n";
        $thisGene->{'Id'} = $sub_cid.".".$thisGene->{'Id'} if !$preserve_ids;
        $thisGene->cropCoords($range->[2], $range->[3]);
        $thisGene->shiftCoords($offset, $scLength, $sub_cid);
        # correcting phase information
        
#        print STDOUT "\t", scalar(@{$gene->{'Exons'}}), " ", scalar(@{$thisGene->{'Exons'}}), "exons left. Coords are shifted by $offset\n";
        if(scalar(@{$thisGene->{'Exons'}}) > 0)  {
#            $gene->dump(\*STDOUT, undef, 1);
            $thisGene->dump(\*STDOUT, undef, 1);
        }
        $index++;
    }
}
