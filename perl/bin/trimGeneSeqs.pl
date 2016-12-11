#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Gene;
use GeneModelProps;
use Organism;
use Contig;
use Genefinder;
use Getopt::Long;
use strict;

my $trimVal = 2000;
my $contigFile = "";
my $outContigs = "";
my $locs = 0;
my $usage = "usage : $0 [-l] -outC outContigs -ctg contigs -d <trimVal> < fastaFile > newFastaFile\nTrim bases from the beginning and end of each contig so that the distance to the beginning or end of a gene would be no more than <trimVal> bases. If -l is specified, outContigs will have the new contigs as locations.\n";

&GetOptions("ctg=s" => \$contigFile,
            "outC=s" => \$outContigs,
            "d=i" => \$trimVal,
            "l" => \$locs
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

use enum qw(STRAND_FWD=0 STRAND_COMP=1);
my $org = Organism->new();

print STDERR "Loading sequences..\n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(OUTCONTIGS, ">$outContigs");

my %contigs = % {$org->loadContigs(\*CONTIG)};
my %genes = % { $org->loadGenes(\*STDIN, 'TRAIN', 1)};
&GenUtils::trimStopCodon(\%genes);

foreach my $cId (keys %contigs) {
    my $contig = $contigs{$cId};
    my $new_cid = $cId;
    my $offset = 1;
    if($cId =~ /^([^\.]+)\.(\d+)$/) {
        $new_cid = $1;
        $offset = $2;
    }
    my $genes = $contig->{'Features'}->{'TRAIN'};
    my @sortedGenes = sort {Gene::compareTo($a, $b);} values %$genes;
    my $min = $contig->{'Len'} + 1;
    my $max = -1;
    my $gene;
    my $exons;
    foreach $gene (@sortedGenes) {
        my $minGene = $gene->lowestCoord();
        my $maxGene = $gene->highestCoord();

        $exons = $gene->{'Exons'};
        if($min > $minGene) {
            $min = $minGene;
        }
        if($max < $maxGene) {
            $max = $maxGene;
        }
    }
    $min = (($min - $trimVal) < 1 ? 1 :  $min - $trimVal);
    $max = (($max + $trimVal) > $contig->{'Len'} ? $contig->{'Len'} : $max + $trimVal);
#    print STDERR "****$contig->{'Id'} $min $max $contig->{'Len'}\n";

    $new_cid = $new_cid.".".($offset + $min - 1);
    my $old_cid = $contig->{'Id'};
    $contig->{'Id'} = $new_cid;
    $contigs{$new_cid} = $contig;

    if($locs) {
        print OUTCONTIGS ">$new_cid\t$old_cid",
        "_", $min,"_", $max, "\n";
    }
    else {
        my $newContigSeq = substr($contig->{'Seq'}, $min - 1, $max - $min + 1);
        print OUTCONTIGS ">$new_cid\n$newContigSeq\n";
    }
    foreach $gene (@sortedGenes) {
        $gene->shiftCoords(-$min + 1, $max - $min + 1, $contig->{'Id'});
        $gene->dump(\*STDOUT, \%contigs, 1);
    }
}
close(OUTCONTIGS); 
close(CONTIG);
