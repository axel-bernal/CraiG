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
my $outDir = "";
my $GCClasses = "0";
my $usage = "usage :$0 -GC gcClasses -ctg contigs -sc sigScoreFile -dir <outDir> < FastaFile\n";

&GetOptions("sc=s" => \$scoreFile,
	    "ctg=s" => \$contigFile,
	    "GC=s" => \$GCClasses,
	    "dir=s" => \$outDir,
	    );

use enum qw(START=0 STOP=1 DONOR=2 ACEPTOR=3);
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
use enum qw(COD NONCOD);

($scoreFile ne "") || die "score file must be specified. ".$usage;
($outDir ne "") || die "output directory must be specified. ".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;
if(!-e $outDir) {
    `mkdir $outDir`;
}
my $i;
my $model = GeneModelProps->new();
my $org = Organism->new();
my $GF = Genefinder->new();
print STDERR "Loading sequences..\n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
my %genes = % { $org->loadGenes(\*STDIN)};
close(CONTIG);
my %gcClasses = %{$model->setGCClasses($GCClasses)};

print STDERR "Done!\nLoading scores...\n";
$GF->retrieveContigRegScores($scoreFile, \%gcClasses, $org->{'Contigs'});
$GF->retrieveContigSigScores($scoreFile,$org->{'Contigs'});
#print join("\n", %{$GF->{'ContigSigScores'}->{'HUMTRPY1B'}}),"\n";
print STDERR "Done!\nPrinting Table\n";

foreach my $id (keys %genes) {
    open(TMP, ">$outDir/$id") || die "$id\n";
    print "$id\n";
    my @score;
    my $currLen = 0;
    my $gene = $genes{$id};
    my $contig = $gene->{'Exons'}->[0]->[0];
    my $strand = $gene->{'Strand'};
    my $i = 0;
    my $p = 0;
    my $exon;
    if($strand == STRAND_FWD) {
	$exon = $gene->{'Exons'}->[$i];
    }
    else {
	$i = scalar(@{$gene->{'Exons'}}) - 1;
	$exon = $gene->{'Exons'}->[$i];
    }
    my @printArray = (0,0,0,0,0,0,0,0);
    for($p = 1; $p <= $contigs{$contig}->{'Len'}; $p++) {
	for(my $j = START; $j <= ACEPTOR; $j++) {
	    $printArray[$j] = 0;
	}
	my $fixed_p = $p + $contigs{$contig}->{'cBeg'} - 1;
	$printArray[4 + $fixed_p%3] = $GF->{'ContigRegScores'}->{$contig}->{$fixed_p}->[$strand]->[COD];
	if(defined($GF->{'ContigSigScores'}->{$contig}->{$fixed_p}->[$strand]->[0])) {
	    my @content = @{$GF->{'ContigSigScores'}->{$contig}->{$fixed_p}->[$strand]};
	    $printArray[$content[0]] = $content[1];
	}
	if($strand == STRAND_FWD) {
	    if($p >= $exon->[1] && $p <= $exon->[2]) {
		$printArray[7] = 1;
	    }
	    elsif($p > $exon->[2]) {
		$printArray[7] = 0;
		if($i < scalar(@{$gene->{'Exons'}}) - 1) {
		    $exon = $gene->{'Exons'}->[++$i];
		}
	    }
	}
	else {
	    if($p >= $exon->[2] && $p <= $exon->[1]) {
		$printArray[7] = 1;
	    }
	    elsif($p > $exon->[1]) {
		$printArray[7] = 0;
		if($i > 0) {
		    $exon = $gene->{'Exons'}->[--$i];
		}
	    }	    
	}
	if(!($p%3) || defined($GF->{'ContigSigScores'}->{$contig}->{$fixed_p}->[$strand]->[0])) {
	    print TMP $p," ", join(" ", @printArray),"\n";    
	}
    }
    close(TMP);
}
print STDERR "Done!\n";
