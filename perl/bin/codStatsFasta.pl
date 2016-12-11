#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Genefinder;
use Contig;
use Getopt::Long;
use strict;

my $scoreFile = "";
my $contigFile = "";
my $GCClasses = "0";
my $usage = "usage : $0 -GC gcClasses -ctg contigs -sc scoreFile < fastaFile\n";

&GetOptions("sc=s" => \$scoreFile,
	    "ctg=s" => \$contigFile,
	    "GC=s" => \$GCClasses,
	    );


($scoreFile ne "") || die "score file must be specified. ".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

use enum qw(EXON=0 INTRON=1 INTERGENIC=4);

my $model = GeneModelProps->new();
my $org = Organism->new();
my $GF = Genefinder->new();
print STDERR "Loading sequences..\n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};
my %gcClasses = %{$model->setGCClasses($GCClasses)};

print STDERR "Done!\nLoading scores...\n";
$GF->retrieveContigRegScores($scoreFile, \%gcClasses, $org->{'Contigs'});
print STDERR "Done!\nPrinting Table\n";
my %locsByGCClass = % {$org->getGCClasses('TRAIN', \%gcClasses)};

foreach my $class (keys %locsByGCClass) {
    my $id;
    foreach $id (@ {$locsByGCClass{$class}}) {
	my $gene = $genes{$id};
	my @scores;
	my $cScorePhase;
	my $codDiffCoding = 0;
	my $codDiffNonCoding = 0;
	my $currLen = 0;
	my $exon = undef;
	my $lastCodPos = 1;
	for(my $i = 0; $i < scalar(@{$genes{$id}->{'Exons'}}); $i++) {
	    my $frame = $currLen % 3;
	    $exon = $gene->{'Exons'}->[$i];
	    my $nonCodType = INTRON;
	    if($i == 0) {
		$nonCodType = INTERGENIC;
	    }
	    my $nonCodDScore = $GF->scoreContigReg($gene->{'Contig'}, $lastCodPos, $exon->[1] - 1, $nonCodType, 0);
	    my $codDScore = $GF->scoreContigReg($gene->{'Contig'}, $exon->[1], $exon->[2], $nonCodType, (3 - $frame)%3);
	    $cScorePhase = $GF->scoreContigReg($gene->{'Contig'}, $exon->[1], $exon->[2], EXON, (3 - $frame)%3);
	    my $cScore = 0;
	    for(my $j = 0;  $j < 3; $j++) {
		$cScore += $GF->scoreContigReg($gene->{'Contig'}, $lastCodPos, $exon->[1] - 1, EXON, $j);
	    }
	    
	    $codDiffCoding += ($cScorePhase - $codDScore)/($exon->[2] - $exon->[1] + 1);
	    if($exon->[1] != $lastCodPos) {
		$codDiffNonCoding += ($cScore/3.0 - $nonCodDScore)/($exon->[1] - $lastCodPos);
	    }	    
	    $currLen = $currLen + $exon->[2] - $exon->[1] + 1;
	    $lastCodPos = $exon->[2] + 1;
	}
	print " " , $id;
	print " " , $class;
	print " ", (sprintf "%.5f", $codDiffCoding);
	print " ", (sprintf "%.5f", $codDiffNonCoding);
	print " ", (sprintf "%.5f", $codDiffCoding - $codDiffNonCoding);
	print "\n";
    }
}
print STDERR "Done!\n";

