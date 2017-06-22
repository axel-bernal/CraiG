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
my $usage = "usage :$0 -GC gcClasses -ctg contigs -sc sigScoreFile -prefix prefix_files\n";
my $prefix = "";

&GetOptions("sc=s" => \$scoreFile,
	    "ctg=s" => \$contigFile,
	    "GC=s" => \$GCClasses,
	    "prefix=s" => \$prefix,
	    );

use enum qw(COD=0 NONCOD=1);
use enum qw(START=0 STOP=1 DONOR=2 ACEPTOR=3);
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

($scoreFile ne "") || die "score file must be specified. ".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $i;
my $window = 54;
my @filenames = ($prefix."starts", $prefix."stops", $prefix."donors", $prefix."aceptors");

my $model = GeneModelProps->new();
my $org = Organism->new();
my $GF = Genefinder->new();
print STDERR "Loading sequences..\n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %gcClasses = %{$model->setGCClasses($GCClasses)};

print STDERR "Done!\nLoading coding scores...\n";
$GF->retrieveContigRegScores($scoreFile, \%gcClasses, $org->{'Contigs'});
print STDERR "Done!\nLoading signal scores...\n";
$GF->retrieveContigSigScores($scoreFile,$org->{'Contigs'});

print STDERR "Done!\nPrinting Table\n";

for(my $sigType = START; $sigType <= ACEPTOR; $sigType++) {
#for(my $sigType = DONOR; $sigType <= DONOR; $sigType++) {
    my %trainingSignals = ();
    $GF->filterSignalsForTraining($org, 1, $sigType, \%trainingSignals, undef, 1);
    open(TMP1, ">$filenames[$sigType]"."-true") || die "$filenames[$sigType]"."-true\n";
    open(TMP2, ">$filenames[$sigType]"."-fals") || die "$filenames[$sigType]"."-fals\n";
    my ($typeDs, $typeUps) = (NONCOD, COD);

    if($sigType == START || $sigType == ACEPTOR) {
	($typeDs, $typeUps) = (COD, NONCOD);
    }
    foreach my $sigL (keys %trainingSignals) {
	my $tmpFile = \*TMP2;
	my $signal = $trainingSignals{$sigL};
	if($signal->{'Strand'} == STRAND_COMP) {
	    next;
	}
	$sigL =~ /(\S+)_(\d+)/ || die "$sigL\n";
	my @score;
	my $cKey = $1;
	my $pos = $2;
	my (@codingScore, @nonCodingScore) = ((0, 0), (0, 0));
	my $maxInd = 0;
	
        for(my $j = 0;  $j < 3; $j++) {
	    if($typeDs == COD) {
		$score[$j] = $GF->scoreContigReg($org, $cKey, $pos + $j + 3, $pos + $window + $j, $typeDs);
	    }
	    else {
		$score[$j] = $GF->scoreContigReg($org, $cKey, $pos + $j - $window, $pos + $j - 3, $typeUps);
	    }
	    if($score[$j] > $score[$maxInd]) {
		$maxInd = $j;
	    }
	}
	$codingScore[COD] = $score[$maxInd];
	if($typeDs == COD) {
	    $codingScore[NONCOD] = $GF->scoreContigReg($org, $cKey, $pos + $maxInd - $window, $pos + $maxInd - 3, $typeDs);
	    $nonCodingScore[COD] = $GF->scoreContigReg($org, $cKey, $pos + 3, $pos + $window, $typeUps);
	    $nonCodingScore[NONCOD] = $GF->scoreContigReg($org, $cKey, $pos - $window, $pos - 3,  $typeUps);
	}
	else {
	    $codingScore[NONCOD] = $GF->scoreContigReg($org, $cKey, $pos + $maxInd + 3, $pos + $window + $maxInd, $typeUps);
	    $nonCodingScore[COD] = $GF->scoreContigReg($org, $cKey, $pos - $window, $pos - 3, $typeDs);
	    $nonCodingScore[NONCOD] = $GF->scoreContigReg($org, $cKey, $pos + 3, $pos + $window, $typeDs);
	}
	
	if($signal->{'Class'} == 1) {
	    $tmpFile = \*TMP1;
	}

	my $gcp = SigUtils::ocurrences($contigs{$cKey}->getSequence(0,$pos - $window, $pos + $window, 0, 0, $signal->{'Strand'}), 'GC',2,1)*1.0/(2*$window);
	my $gcc = SigUtils::ocurrences($contigs{$cKey}->getSequence(0,$pos - $window, $pos + $window, 0, 0, $signal->{'Strand'}), 'G|C',1,1)*1.0/$window;
	
	print $tmpFile ">$cKey"."_".$signal->{'Strand'}."_"."$pos ", (sprintf "%.1f", $codingScore[COD]);
	print $tmpFile " ", (sprintf "%.1f", $codingScore[NONCOD]);
	print $tmpFile " ", (sprintf "%.1f", $nonCodingScore[COD]);
	print $tmpFile " ", (sprintf "%.1f", $nonCodingScore[NONCOD]);
        print $tmpFile " ", (sprintf "%.5f", $gcp);
        print $tmpFile " ", (sprintf "%.5f", $gcc);
	print $tmpFile " ", $signal->{'Score'};
	print $tmpFile "\n";
    }
    close(TMP1);close(TMP2);
    `cat "$filenames[$sigType]-true" "$filenames[$sigType]-fals" > $filenames[$sigType]`;
}

print STDERR "Done!\n";
