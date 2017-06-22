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
my $usage = "usage :$0 -GC gcClasses -ctg contigs -sc sigScoreFile < fastaFile\n";

&GetOptions("sc=s" => \$scoreFile,
	    "ctg=s" => \$contigFile,
	    "GC=s" => \$GCClasses,
	    );

use enum qw(COD=0 NONCOD=1);
use enum qw(START=0 STOP=1 DONOR=2 ACEPTOR=3);
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

($scoreFile ne "") || die "score file must be specified. ".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

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

open(TMP1,">starts"); open(TMP2,">stops");open(TMP3,">donors");open(TMP4,">aceptors");

foreach my $id (keys %genes) {
    print ">$id $genes{$id}->{'Strand'}\n";
    my @score;
    my $currLen = 0;
    my $gene = $genes{$id};
    for($i = 0; $i < scalar(@ {$gene->{'Exons'}}); $i++) {
	my $exon = $gene->{'Exons'}->[$i];
        my $frame  = $currLen % 3;
	my ($posSig1, $posSig2)= ($exon->[1], $exon->[2]);
	my ($sig1, $sig2) = (ACEPTOR, DONOR);
	my ($shift1, $shift2) = $gene->{'Strand'} == STRAND_FWD?(-2, 1):(2,-1);
	my ($tmp1, $tmp2) = (\*TMP4, \*TMP3);
        for(my $j = 0;  $j < 3; $j++) {
            $score[($j + $frame)%3] = $GF->scoreContigReg($exon->[0], $exon->[1] + $j, $exon->[2] + $j, COD);
        }
        $currLen = $currLen + $exon->[2] - $exon->[1] + 1;
	if($i == 0) {
	    $sig1 = START; $tmp1 = \*TMP1;
	}
	else {
	    $posSig1 += $shift1;
	}
	if($i == scalar(@ {$gene->{'Exons'}}) - 1) {
	    $sig2 = STOP; $tmp2 = \*TMP2;
	}
	$posSig2 += $shift2;

	my $content1 = $GF->{'ContigSigScores'}->{$exon->[0]}->{$posSig1}->[$gene->{'Strand'}];
	my $content2 = $GF->{'ContigSigScores'}->{$exon->[0]}->{$posSig2}->[$gene->{'Strand'}];
	($content1->[0] == $sig1 && $content2->[0] == $sig2 && $content1->[2] == $content2->[2] && $content2->[2] == 1) || die "$posSig1 $content1->[2] $sig1 $posSig2 $content2->[2] $sig2\n";

	my $gcp = SigUtils::ocurrences($contigs{$exon->[0]}->getSequence(0,$exon->[1], $exon->[2], 0, 0, $gene->{'Strand'}), 'GC',2,1)*1.0/(abs($exon->[1] - $exon->[2]) + 1);
	my $gcc = SigUtils::ocurrences($contigs{$exon->[0]}->getSequence(0,$exon->[1], $exon->[2], 0, 0, $gene->{'Strand'}), 'G|C',1,1)*1.0/(abs($exon->[1] - $exon->[2]) + 1);

	print $tmp1 " ", $content1->[1];
        print $tmp1 " ", (sprintf "%.5f", $gcc);
        print $tmp1 " ", (sprintf "%.5f", $gcp);
        print $tmp1 " ", abs($exon->[1] - $exon->[2]) + 1, " $frame";
	for(my $j = 0;  $j < 3; $j++) {
            print $tmp1 " ", (sprintf "%.5f", $score[$j]);
        }
	print $tmp1 "\n";
	
	print $tmp2 " ", $content2->[1];
        print $tmp2 " ", (sprintf "%.5f", $gcc);
        print $tmp2 " ", (sprintf "%.5f", $gcp);
        print $tmp2 " ", abs($exon->[1] - $exon->[2]) + 1, " $frame";
	for(my $j = 0;  $j < 3; $j++) {
            print $tmp2 " ", (sprintf "%.5f", $score[$j]);
        }
	print $tmp2 "\n";
    }
}
close(TMP1);close(TMP2);close(TMP3);close(TMP4);
print STDERR "Done!\n";
