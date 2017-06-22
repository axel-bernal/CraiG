package SigUtils;
require Exporter;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
use Utils;
@ISA = ( Exporter );
@EXPORT= qw (
	     getSeqClassBySignal
	     getGClass
	     getContigClustersBySignal
	     ocurrences
	     regIndex
	     signalsToSVMStartFormat
	     signalsToSVMDonorFormat
	     signalsToSVMUnivFormat
	     signalToSVMUnivFormat
	     testSignals
	     storeSignals
	     trainSignals
	     );

use IO::Handle;
use strict;

### Enumerates ###

use enum qw(START STOP DONOR ACEPTOR);
use enum qw(GC_CLUSTER_BEG=0 GC_CLUSTER_END=1 GC_CLUSTER_CLASS=2 GC_CLUSTER_SEQ=3);
use enum qw(COMP=4);
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my %condEncoding  = ('AA' => 0,
		  'AC' => 1,
		  'AG' => 2,
		  'AT' => 3,
		  'CA' => 4,
		  'CC' => 5,
		  'CG' => 6,
		  'CT' => 7,
		  'GA' => 8,
		  'GC' => 9,
		  'GG' => 10,
		  'GT' => 11,
		  'TA' => 12,
		  'TC' => 13,
		  'TG' => 14,
		  'TT' => 15,
		  );    

my %uniEncoding = ('A' => [0, 0, 0, 1],
		'C' => [0, 0, 1, 0],
		'T' => [0, 1, 0, 0],
		'G' => [1, 0, 0, 0],
		);

sub getContigClustersBySignal {
    my($contigs, $sigClasses, $window, $signal, $lenSignal, $step) = @_;
    my %contigClusters = ();
    defined($step) or $step = $window;

    foreach my $cKey (keys %$contigs) {
	my $contig = $contigs->{$cKey};
	my $base = $contig->{'Beg'};
	my $init = $base;
	my $seqLen;
	my $lastInsComp = undef;
	my $lastInsFwd = undef;
	while($base <= $contig->{'End'}) {
	    my $currWindow = $window;
	    my $seq = $contig->getSequence(0, $base, $base + $currWindow - 1, 0, 1, STRAND_FWD);
	    my $compSeq = &GenUtils::revComp($seq);
	    my $seqLen = length($seq);
	    my $seqClass = &SigUtils::getSeqClassBySignal($signal, $lenSignal, $seq, $sigClasses);
	    my $compSeqClass = &SigUtils::getSeqClassBySignal($signal, $lenSignal, $compSeq, $sigClasses);
	    
#	    print "$base $seqClass $compSeqClass ",length($seq), "\n";
	    if(defined($lastInsFwd) && ($lastInsFwd->[GC_CLUSTER_CLASS] == $seqClass)) {
		$lastInsFwd->[GC_CLUSTER_END] += $step;
#		$lastInsFwd->[GC_CLUSTER_SEQ] = $lastInsFwd->[GC_CLUSTER_SEQ].$seq;
#		print "increased $lastInsFwd->[GC_CLUSTER_BEG] $lastInsFwd->[GC_CLUSTER_END] $lastInsFwd->[GC_CLUSTER_CLASS] \n";
	    }
	    else {
#		print "added $base ", $base + $currWindow - 1, " $seqClass ", length($seq), "\n";
		push(@{$contigClusters{$contig->{'Id'}}->[STRAND_FWD]}, [$base, $base + $seqLen - 1, $seqClass]);
		$lastInsFwd = $contigClusters{$contig->{'Id'}}->[STRAND_FWD]->[scalar(@{$contigClusters{$contig->{'Id'}}->[STRAND_FWD]}) - 1];
	    }
	    
	    if(defined($lastInsComp) && ($lastInsComp->[GC_CLUSTER_CLASS] == $compSeqClass)) {
		$lastInsComp->[GC_CLUSTER_BEG] += $step;
#		$lastInsComp->[GC_CLUSTER_SEQ] = $compSeq.$lastInsComp->[GC_CLUSTER_SEQ];
#		print "increased $lastInsComp->[GC_CLUSTER_BEG] $lastInsComp->[GC_CLUSTER_END] $lastInsComp->[GC_CLUSTER_CLASS] \n";
	    }
	    else {
#		print "added ",2*$contig->{'Beg'} + $contig->{'Len'} - 1 - $base, " ", 2*$contig->{'Beg'} + $contig->{'Len'} - ($base + $currWindow), " $seqClass ", length($compSeq), "\n";
		unshift(@{$contigClusters{$contig->{'Id'}}->[STRAND_COMP]}, [$base + $seqLen - 1, $base, $compSeqClass]);
		$lastInsComp = $contigClusters{$contig->{'Id'}}->[STRAND_COMP]->[0];
	    }
	    $base += $step;
	}
    }
    return \%contigClusters
}

sub getGCPerc {
    my($seq) = @_;
# Ad-hoc code for improving performance
    my $countGC = 0.0;
    for(my $i = 0; $i < length($seq); $i++) {
	if(substr($seq, $i, 1) =~ /g/i || substr($seq, $i, 1) =~ /c/i) {
	    $countGC++;
	}
    }
    return $countGC*100.0/length($seq);
}

sub getWindowGCPerc {
    my($seq, $window) = @_;
# Ad-hoc code for improving performance
    my $ctxc = 0.0;
    my $i;
    my @cts;
    if(!$window || length($seq) < $window) {
        $window = length($seq);
    }

    # 5 prime

    for($i = 0; $i < $window; $i++) {
	if(substr($seq, $i, 1) =~ /g/i || substr($seq, $i, 1) =~ /c/i) {
	    $ctxc++;
	}
    }
    
    # middle
    for($i = $window/2; $i < length($seq) - $window/2; $i++) {
	if($i % ($window/2)) {
	    push(@cts, ($ctxc*100/$window));
	}

	if(substr($seq, $i - $window/2, 1) =~ /g/i || 
	   substr($seq, $i - $window/2, 1) =~ /c/i) {
	    $ctxc--;
	}
	
	if(substr($seq, $i + $window/2, 1) =~ /g/i || 
	   substr($seq, $i + $window/2, 1) =~ /c/i) {
	    $ctxc++;
	}

    }
    
    return \@cts;
}


sub getGCClass {
    my($seq, $gcClasses) = @_;
# Ad-hoc code for improving performance
    my $idClass = 100;
    my $sigContent = getGCPerc($seq);

    foreach my $class (keys %$gcClasses) {
	my $range = $gcClasses->{$class};
	if($sigContent >= $range->[0] && $sigContent <= $range->[1]) {
	    $idClass = $class;
	    last;
	}
    }

    return $idClass;
}

sub getSeqClassBySignal {
    my ($signal, $lenSignal, $seq, $sigClasses) = @_;
    my $class;
    my $sigContent = ocurrences($seq, $signal, $lenSignal)*100.0/length($seq);

#    print "GC $signal, $lenSignal, ", ocurrences($seq, $signal, $lenSignal), " ", length($seq)," $sigContent\n";
    # getting the respective GC content class
    foreach $class (keys %$sigClasses) {
	my $range = $sigClasses->{$class};
	if($sigContent >= $range->[0] && $sigContent <= $range->[1]) {
	    return $class;
	}
    }
    return "";
}

# Get ocurrences of pattern in seq *
sub ocurrences {
    my($seq, $pattern, $lenPattern, $step) = @_;
    my($base,$numMatches) = (0,0);
#    print $seq,"\t",$pattern,"\t";
    if(!defined($lenPattern)) {
	$lenPattern = length($pattern);
    }
    $step = 1 if ! defined($step);
    my $initPos = regIndex($seq, $pattern, $lenPattern, $base, $step);
    while($initPos >= 0) {
        $numMatches++;
        $base = $initPos;
        $initPos = regIndex($seq, $pattern, $lenPattern, $base + 1, $step);
    }
    return $numMatches;
}

sub regIndex {
    my($seqdna, $expr, $len, $base, $step) = @_;
    my($tmpstr) = "";
    while(1) {
        if($base + $len > length($seqdna)) {
            return -1;
        }
        $tmpstr = substr($seqdna, $base, $len);
        if($tmpstr =~ /^($expr)/i) {
            return $base;
        }
        else {
            $base += $step;
        }
    }
}

sub getBinaryRepSeqVector {
    my($seq) = @_;
    my @tmp = ();
    my @dummy = (0,0,0,0);
    for(my $i = 0; $i < length($seq); $i++) {
        my $character = substr($seq,$i,1);
#       print $character, " ", $$encoding{$character},"\n";
        if(!defined($uniEncoding{$character})) {
            push(@tmp, @dummy);
        }
        else {
            push(@tmp , @ { $uniEncoding{$character}});
        }
    }
    return \@tmp;
}

sub signalsToSVMStartFormat {
    my($signals, $patterns, $outFile) = @_;
    open(TMP,">$outFile") || die "can't open $outFile\n";;
    my $signal;
    my $i;
    foreach $signal (@$signals) {
	my $numMatches;
	my @vector = ();
	my $pattern;
#	print join(" ",%$signal),"\n";
	foreach $pattern ( @$patterns) {
	    $numMatches = &ocurrences($signal->{'Seq'},$pattern);
#	    print "PAT ", join(" ", @$pattern)," $numMatches\n";
	    push(@vector,$numMatches);
	}
	print TMP "$signal->{'Class'}";
	for($i=0;$i < scalar(@vector);$i++) {
	    print TMP " ", $i+1,":",$vector[$i];
#       print TMP " ", $vector[$i];
	}
	print TMP "\n";
    }
    close(TMP);
    return;
}


sub signalsToSVMDonorFormat {
    my ($signals, $fileHandle) = @_;
    open(TMP,">$fileHandle") || die "can't open $fileHandle\n";;
    my $signal;
    my $ind;
    my $j;
    my $i;
    my $base;
    my @vector = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    foreach $signal (@$signals) {
	print TMP "$signal->{'Class'}";
	$base = 0;
	for($i = 0; $i < length($signal->{'Seq'}) - 1; $i++) {
	    $ind = $condEncoding{substr($signal->{'Seq'}, $i, 2)};
#	    print substr($signal->{'Seq'}, $i, 2)," ",$ind,"\t";
	    if(defined($ind)) { $vector[$ind] = 1; }
	    for($j=0;$j < scalar(@vector); $j++) {
		print TMP " ", $j+$base+1,":",$vector[$j];
	    }
	    if(defined($ind)) { $vector[$ind] = 0; }
	    $base += $j;
	}
	print TMP "\n";
    }
    close(TMP);
    return;
}

sub signalToSVMUnivFormat {
    my($signal, $outFile) = @_;
    my @vector = @{&getBinaryRepSeqVector($signal->{'Seq'})};
    print $outFile "$signal->{'Class'}";
    for(my $i = 0; $i < scalar(@vector); $i++) {
	print $outFile " ",$i+1,":",$vector[$i];
    }
    print $outFile "\n";
}

sub signalsToSVMUnivFormat {
    my($signals, $outFile) = @_;
    my $signal;
    open(TMP,">$outFile") || die "can't open $outFile\n";;
    TMP->autoflush(1);
    foreach $signal (@$signals) {
	signalToSVMUnivFormat($signal, \*TMP);
    }
    close(TMP);
    return;
}

sub testSignals {
    my($signals, $modelFile, $signalFile, $sigType, $tag, $GF) = @_;
    my $predictions = "$signalFile"."$tag.test.pred";    
    my $i = 0;
    my @results = ();
    my $motifFlags = " -w 8 -l 2 ";
    my @patterns = ();
    my @sigNames = ("START","STOP","DONOR","ACEPTOR");
    storeSignals($signals, "$signalFile"."$tag.test.fsa");
    signalsToSVMUnivFormat($signals, $signalFile."$tag.test.svm");
    if($sigType == START) {
	system("score_contigs_nn -t $sigNames[$sigType]"."_NET -f $signalFile"."$tag.test.fsa | egrep -v \"0\$\" | grep -v \"\>\" > $predictions");
    }
    else {
	system("svm_classify $signalFile"."$tag.test.svm $modelFile $predictions");
    }


#    if($sigType == DONOR) {
#	signalsToSVMDonorFormat($signals, $signalFile."$tag.test.svm");
#    }
#    elsif($sigType == START) {
#	system("motif $motifFlags -k ".(scalar(@$signals)*0.65)." -f $signalFile"."$tag.test.fsa > $signalFile.test.pat");
#	@patterns = @ {&Utils::fileToArray($signalFile.".test.pat")};
#	signalsToSVMStartFormat($signals, \@patterns, $signalFile."$tag.test.svm");
#	unlink("$signalFile.test.pat");
#    }
#    else {
#    signalsToSVMUnivFormat($signals, $signalFile."$tag.test.svm"); 	
#    }

#    system("svm_classify $signalFile"."$tag.test.svm $modelFile $predictions");
    @results = @ {&Utils::fileToArray($predictions)};
#    unlink $predictions;
    return \@results;
}

sub storeSignals { 
    my ($signals, $file, $cat) = @_;
    my  $signal;
    if(defined($cat) && $cat == 1) {
	open(TMP,">>$file") || die "Can't open $file\n";
    }
    else {
	open(TMP,">$file") || die "Can't open $file\n";
    }
    foreach $signal (@$signals) {
	$signal->{'Pos'} =~ /(\S+)_(\d+)$/;
        print TMP ">$1"."_".$signal->{'Strand'}."_"."$2\t$signal->{'Class'}\n",$signal->{'Seq'},"\n";
    }
    close(TMP);
}

sub trainSignals {
    my($signals, $modelFile, $signalFile, $learnFlags, $sigType, $GF) = @_;
    my $signal;
    my @patterns;
    my $motifFlags = " -w 8 -l 2 ";
    storeSignals($signals, "$signalFile.fsa");
#    if($sigType == DONOR) {
#	signalsToSVMDonorFormat($signals, "$signalFile.svm"); 
#    }
#    elsif($sigType == START) {
#	system("motif $motifFlags -k ".(scalar(@$signals)*0.65)." -f $signalFile.fsa > $signalFile.pat");
#	@patterns = @ {&Utils::fileToArray($signalFile.".pat")};
#	signalsToSVMStartFormat($signals, \@patterns, "$signalFile.svm");
#    }
#    else {
	signalsToSVMUnivFormat($signals, "$signalFile.svm"); 	
#    }
    system("svm_learn $learnFlags $signalFile.svm $modelFile");
}
1;

