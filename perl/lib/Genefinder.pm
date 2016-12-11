package Genefinder;
require Exporter;
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Histogram;
use Contig;
use Signal;
use Gene;
use GeneModelProps;
@ISA = ( Exporter );
@EXPORT= qw (
	     viterbiContig
	     viterbiGene
	     computeIntergenicHistogram
	     computeHistograms
	     extractGeneSignals
	     addNewTrainingSignals
	     trainContigSignals
	     trainCodingRegions
	     trainIStateRegions
	     trainGHMMParams
	     scoreContigSignals
	     genGCIsochores
	     scoreGCIsochores
	     genContigClusters
	     scoreContigClusters
	     setContigSignals
	     classifySignalsWGeneSet
	     filterSignalsForTraining
	     retrieveContigSigScores
	     storeContigSigScores
	     retrieveContigRegScores
	     scoreContigReg
	     storeClusters
	     );

use strict;
### Enumerates ###
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
use enum qw(GC_CLUSTER_BEG=0 GC_CLUSTER_END=1 GC_CLUSTER_CLASS=2 GC_CLUSTER_SEQ=3);
use enum qw(START=0 STOP=1 DONOR=2 ACEPTOR=3);
use enum qw(COMP=4);
use enum qw(TR_SIG_TYPE=0 TR_PRED_ST=1 TR_CONTIG=2 TR_POS=3 TR_SIG_SCORE=4 TR_PTR_NEXT=5); # transition fields.
use enum qw(ROOT=-2 INIT=-1 EXON=0 INTRON=1 INTERGENIC=4); # states, including the synchronizing ones
use enum qw(BEGS=8 ENDS=9); # additional synchronizing transitions
use enum qw(CODING=5 INTRON_LEN=6 EXON_LEN=7 EXON_NUM=8); #params
use enum qw(PATH_SCORE=0 PATH_NEXT_ST=1 PATH_LAST_ST=2);

sub new {
    my $type = shift;
    my %params = @_;
    
    my $self = {
	'Model' => $params{'Model'},
	'PredState' => [INTERGENIC, EXON, EXON, INTRON,
			INTERGENIC, EXON, INTRON, EXON,
			INIT, INTERGENIC],
	'HistExonLen' => undef,
	'HistExonNum' => undef,
	'HistIntronLen' => undef,
	'HistIntergenicLen' => undef,
	'Params' => undef,
	'ContigGCClusters' => undef,
	'ContigRegScores' => undef,
	'ContigSigScores' => undef, # array of signals STARTS, STOPS, DONORS, ACEPTORS and complementary signals, existent in the contigs
    };
    bless $self, $type;
    return $self;
}

#### Recursive method that will insert each new transition in our tree of possible valid parses of the sequence
sub insertTransition {
    my($trans, $entry, $minExonLen, $maxExonLen, $minIntronLen, $maxIntronLen, $maxIntergenicLen, $maxNumExons, $currNumExons, $currCodingLen) = @_;
    my $nextTrans;
    my $insertIt;
    my $insertedNext;
    my $alreadyInserted = 0;
    my $numAddExons = 0;
    my $addCodingLen = 0;
#    print "Matching $entry->[TR_SIG_TYPE] \@ $entry->[TR_POS] with $trans->[TR_SIG_TYPE] \@ $trans->[TR_POS]\n";
    # Checking if entry has been already inserted at this level
    foreach $nextTrans (@{$trans->[TR_PTR_NEXT]}) {
	if($entry == $nextTrans) {
	    $alreadyInserted = 1;
	    last;
	}
    }

    if(!defined($currNumExons)) {
	$currNumExons = 0;
    }
    $insertIt = 0;
    if(!$alreadyInserted) {
	############ inserting non-coding regions #############
	#after an exon, we have to make sure to insert an INTRON
	if($entry->[TR_PRED_ST] == INTRON && ($entry->[TR_POS] - $trans->[TR_POS]) <= $maxIntronLen && ($entry->[TR_POS] - $trans->[TR_POS]) >= $minIntronLen && ($trans->[TR_SIG_TYPE] == DONOR && $entry->[TR_SIG_TYPE] == ACEPTOR)) {
	    $insertIt = 1;
	}
	#if we receive an ENDS transition, a STOP must be right before it
	elsif($entry->[TR_PRED_ST] == INTERGENIC && ($entry->[TR_POS] - $trans->[TR_POS]) <= $maxIntergenicLen && $trans->[TR_POS] == STOP) {
	    $insertIt = 1;
	}
	############ inserting coding regions ############
	# after an intron we must insert an exon making sure it doesn't break the length distribution and is in frame with oncomming transition
	elsif($entry->[TR_PRED_ST] == EXON && ($entry->[TR_POS] - $trans->[TR_POS]) >= $minExonLen && ($entry->[TR_POS] - $trans->[TR_POS]) <= $maxExonLen && $currNumExons + 1 < $maxNumExons) {
	    # beginning in START ending in a DONOR
	    if($trans->[TR_SIG_TYPE] == START && $entry->[TR_SIG_TYPE] == DONOR) {
		$insertIt = 1;
	    }
	    # beginning in an ACEPTOR ending in an DONOR
	    elsif($trans->[TR_SIG_TYPE] == ACEPTOR && $entry->[TR_SIG_TYPE] == DONOR) {
		$insertIt = 1;
	    }
	    # beginning in an START ending in a STOP
	    elsif($trans->[TR_SIG_TYPE] == START && $entry->[TR_SIG_TYPE] == STOP && $entry->[TR_POS] % 3 == $trans->[TR_POS] % 3) {
		$insertIt = 1;
	    }
	    # beginning in an ACEPTOR ending in an STOP
	    elsif($trans->[TR_SIG_TYPE] == ACEPTOR && $entry->[TR_SIG_TYPE] == STOP) {
		if(!(($currCodingLen + $entry->[TR_POS] - $trans->[TR_POS]) % 3)) {
		    $insertIt = 1;
		}
	    }
	}
	if($insertIt) {
#	    print ">" x $currNumExons, "\@ $currCodingLen and $entry->[TR_SIG_TYPE] @ $entry->[TR_POS] with $trans->[TR_SIG_TYPE] @ $trans->[TR_POS]\n";
	    push(@{$trans->[TR_PTR_NEXT]}, $entry);	
	}
    }
    if(!defined($trans->[TR_PTR_NEXT])) { # empty list of pointers to valid transitions
	my(@emptyList) = ();
	$trans->[TR_PTR_NEXT] = \@emptyList;	    
    }
    
    foreach $nextTrans (@{$trans->[TR_PTR_NEXT]}) {
	$addCodingLen = 0;
	$numAddExons = 0;
	
	### Updating the number of exons we have traversed up to this point
	#### check we have not inserted entry already in this path
	if($nextTrans eq $entry) {
	    next;
	}
	if($nextTrans->[TR_PRED_ST] == EXON) {
	    $numAddExons = 1;
	    $addCodingLen = $nextTrans->[TR_POS] - $trans->[TR_POS];
	    if($trans->[TR_SIG_TYPE] == ACEPTOR) {
		$addCodingLen -= 2;
	    }
	    if($nextTrans->[TR_SIG_TYPE] == STOP) {
		$addCodingLen = -$currCodingLen;
	    }
	}
	if(insertTransition($nextTrans, $entry, $minExonLen, $maxExonLen, $minIntronLen, $maxIntronLen, $maxIntergenicLen, $maxNumExons, $currNumExons + $numAddExons, $currCodingLen + $addCodingLen)) {
	    $insertedNext = 1;
	}
    }
    return $insertIt;
}

### Viterbi algorithm for contigs, it computes the bestN possible paths whenever a super-state (gene) finishes. In this way at the end we will obtain the best parse taking into account the bestN parses up to a certain point in the sequence ######
sub viterbiContig {
    my($self, $contig, $from, $to, $maxIntergenicLen, $minSepSignals, $bestN, $minCoverage) = @_;
    my($i, $counter, $j, $k, $transScore) = (0, 0, 0, 0,0);
    my $strand;
    my $contigId = $contig->{'Id'};
    my $model = $self->{'Model'};
    my $base = $from;
    my @transitions = ();
    my @smallestGeneEnd = ();
    my($minExonLen, $maxExonLen, $maxIntronLen, $minIntronLen, $maxExonNum) = (0,0,0,0,0);
    #### build up a root for traverse purposes ####
    my @transRoot = ();
    my @potentialGenes = ();
    for($strand = STRAND_FWD; $strand <= STRAND_COMP; $strand++) {
	my @lastSignal = ();
	print "strand $strand\n";
	#synchronizing beginning of sequence
	@{$transitions[$strand]} = ();
	$transRoot[$strand] = [BEGS, $self->{'PredState'}->[BEGS], "", $from, 0, $transitions[$strand]];
        #### going through the contigs in  strand $strand
	my $cluster = 0;
	my $currentGCCluster = $self->{'ContigGCClusters'}->{$contigId}->[$strand]->[$cluster];
	while($from > &Utils::max($currentGCCluster->[GC_CLUSTER_BEG], $currentGCCluster->[GC_CLUSTER_END])) {
	    $cluster++;
	    $currentGCCluster = $self->{'ContigGCClusters'}->{$contigId}->[$strand]->[$cluster];
	}
	($to < &Utils::max($self->{'ContigGCClusters'}->{$contigId}->[$strand]->[scalar(@{$self->{'ContigGCClusters'}->{$contigId}->[$strand]}) - 1]->[GC_CLUSTER_BEG], $self->{'ContigGCClusters'}->{$contigId}->[$strand]->[scalar(@{$self->{'ContigGCClusters'}->{$contigId}->[$strand]}) - 1]->[GC_CLUSTER_END])) || die "$to out of bonds\n";
	my $GCClass;	
	for($i = $from; $i <= $to; $i++) {
	    if($i == $from || $i > &Utils::max($currentGCCluster->[GC_CLUSTER_BEG], $currentGCCluster->[GC_CLUSTER_END])) {
		# initializing maximum and minimum values for the GC class associated with this cluster
		$GCClass = $currentGCCluster->[GC_CLUSTER_CLASS];
		$minExonLen = $self->{'HistExonLen'}->{'Min'}-10;
		$maxExonLen = $self->{'HistExonLen'}->{'Max'}+10;
		$maxExonNum = $self->{'HistExonNum'}->{'Max'}*2;
		$maxIntronLen = $self->{'HistIntronLen'}->{$GCClass}->{'Max'}+100;
		$minIntronLen = $self->{'HistIntronLen'}->{$GCClass}->{'Min'}-100;
		$cluster++;
		$currentGCCluster = $self->{'ContigGCClusters'}->{$contigId}->[$strand]->[$cluster];
		print STDERR "switching cluster at $i\n$minExonLen $maxExonLen $maxExonNum $maxIntronLen $minIntronLen\n";
	    }
	    
	    if($i - $base > $maxIntergenicLen + $maxExonNum*($maxExonLen + $maxIntronLen)*2 || $i == $to)  {
		print STDERR "time to check genes at $i\n";
		if($i == $to) {
		    ### end of sequence : synching transition ###
		    my($entry) = [ENDS, $self->{'PredState'}->[ENDS], $contigId, $to, 0, undef];
		    insertTransition($transRoot[$strand], $entry, $minExonLen, $maxExonLen, $minIntronLen, $maxIntronLen, $maxIntergenicLen, $maxExonNum, 0, 0);
		}
		### decode tree returning the bestN possible gene parses ###
		my @genes = ();
		my @currentExons = ();
		$self->viterbiGene(\@genes, \@currentExons, 0, $transRoot[$strand], 0, $bestN, $strand);
		### pruning tree up to the end of the smallest gene found so far for both strands ###
		$smallestGeneEnd[$strand] = $i;
		foreach my $gene (@genes) {
		    if($gene->{'Exons'}->[scalar(@{$gene->{'Exons'}}) - 1]->[2] < $smallestGeneEnd[$strand]) {
			$smallestGeneEnd[$strand] = $gene->{'Exons'}->[scalar(@{$gene->{'Exons'}}) - 1]->[2];
		    }
		}
		print "Smallest gene end = ", $smallestGeneEnd[$strand], "\n";
		while(($transRoot[$strand]->[TR_PTR_NEXT]->[0]->[TR_SIG_TYPE] == START) && ($transRoot[$strand]->[TR_PTR_NEXT]->[0]->[TR_POS] < $smallestGeneEnd[$strand]) && ($i - $transRoot[$strand]->[TR_PTR_NEXT]->[0]->[TR_POS] > ($maxExonLen + $maxIntronLen))) {
		    print "Taking out ", $transRoot[$strand]->[TR_PTR_NEXT]->[0]->[TR_POS],"\n";
		    shift(@{ $transRoot[$strand]->[TR_PTR_NEXT]});
		}
		# returning the complement genes to its original location
		if($strand == STRAND_COMP) {
		    foreach my $gene (@genes) {
			print "normal gene ", join(" ",%$gene);
			$gene->toComplement($to + $from);
			print "complemented gene ", join(" ",%$gene);
		    }
		}
		$base = $i; ## updating base
		push(@potentialGenes, @genes); # pushing genes into array
	    }
	    
	    if(!defined($self->{'ContigSigScores'}->{$contigId}->{$i}->[$strand])) { 
		next;
		print "continue ... $contigId\n";
	    }
	    my($type, $score) = @ {$self->{'ContigSigScores'}->{$contigId}->{$i}->[$strand]};
	    my $entry = [$type, $self->{'PredState'}->[$type], $contigId, $i, $score, undef];
	    
	    if(!defined($lastSignal[$type]->[$i % 3]) || (defined($lastSignal[$type]->[$i % 3]) && ($i - $lastSignal[$type]->[$i % 3]) >= $minSepSignals)) {
#		print "inserted ", join("|", @{$entry}),"\n";
		if($type == START) {
		    push(@ {$transitions[$strand]}, $entry);
		    next;
		}
		
		insertTransition($transRoot[$strand], $entry, $minExonLen, $maxExonLen, $minIntronLen, $maxIntronLen, $maxIntergenicLen, $maxExonNum, 0, 0);
	    }
	    $lastSignal[$type]->[$i % 3] = $i;
	}
    }

    ## stuff used in dynamical programming for the contig
    my (@bestParses) = ();
    my (@matScores) = ();
    my (@btMatrix) = ();
    my ($numGenes,$maxForSt,$initMaxScore) = (0,0,0);
    my @sortedPotGenes = sort Gene::compareTo (@potentialGenes);
    my $numParses = scalar(@sortedPotGenes) > 2*$bestN ? 2*$bestN:scalar(@sortedPotGenes);
    my @maxForSt = ();
    #### Computing the best possible arrangement of genes stored in @sortedPotGenes
    ## initializing matrices

    for($j = 0; $j < $numParses; $j++) {
	$bestParses[$numGenes]->[$j] = $j;
	$matScores[$numGenes]->[$j] = $sortedPotGenes[$bestParses[$numGenes]->[$j]]->{'Score'};
	$btMatrix[$numGenes]->[$j] = -1;
        if($matScores[$numGenes]->[$j] > $initMaxScore) {
            $maxForSt[$numGenes] = $j;
            $initMaxScore = $matScores[$numGenes]->[$j];
        }
        $btMatrix[$numGenes]->[$j] = -1;
    }
    $numGenes++;
    while(1) {
	# put the best possible succesors for the current $numParses in bestParses[$numGenes - 1], put results in @possibleSuccGenes
	my @possibleSuccGenes = @ { &GenUtils::findPossibleGeneSuccesors(\@sortedPotGenes, $bestParses[$numGenes - 1], $maxIntergenicLen*2)};
	my @allParses = ();
	
	for($j = 0; $j < scalar(@possibleSuccGenes); $j++) {
	    my $maxScore = 0;
	    $matScores[$numGenes]->[$j] = -1;
	    $btMatrix[$numGenes]->[$j] = -1;
	    for($i = 0; $i < $numParses; $i++) {
		if($bestParses[$numGenes - 1]->[$i] == $possibleSuccGenes[$j]) {
		    next;
		}
		my @gap = @{$sortedPotGenes[$bestParses[$numGenes - 1]->[$i]]->Gap($possibleSuccGenes[$j])};
		my $gapConsGenes = $gap[1] - $gap[0];
		if($gapConsGenes <= 0 || $gapConsGenes > $maxIntergenicLen) {
		    next;
		}
		my $thisScore = $self->{'Params'}->[INTERGENIC]->[$gapConsGenes/$model->{'HistIntergenicBinLen'}] + $sortedPotGenes[$bestParses[$numGenes - 1]->[$i]]->{'Score'};
		if($thisScore > $maxScore) {
		    $allParses[$j] = [$maxScore, $j, $i];
		    $maxScore = $thisScore;
		}
	    }
	}
	#Ordering allParses
	my @allOrderedParses  = sort Utils::byScore @allParses;
	#Getting bestParses, Scores and btMatrixes
	$numParses = (scalar(@possibleSuccGenes) > 2*$bestN) ? 2*$bestN : scalar(@possibleSuccGenes);	
	for($j = 0; $j < $numParses; $j++) {
	    $bestParses[$numGenes]->[$j] = $possibleSuccGenes[$allParses[$j]->[PATH_NEXT_ST]];
	    $btMatrix[$numGenes]->[$j] = $allOrderedParses[$j]->[PATH_LAST_ST];
	    $matScores[$numGenes]->[$j] = $allOrderedParses[$j]->[PATH_SCORE];
	}
	#Getting maximum of last iteration
	my $maxScore = 0;
        $maxForSt[$numGenes] = -1;
        for($j = 0; $j < $numParses; $j++) {
            if($matScores[$numGenes]->[$j] > $maxScore) {
                $maxScore = $matScores[$numGenes]->[$j];
                $maxForSt[$numGenes] = $j;
            }
        }                           
        #### checking stop condition #####
	if($maxForSt[$numGenes] < 0) { # I couldnt find a good continuation, leaving
            last;
        }
        if($btMatrix[$numGenes]->[$maxForSt[$numGenes]] < 0) { # We are at the beginning
            $numGenes++;
            next;
	}
	$numGenes++;
    }

    #Check each of the maximum values reached at each numGene
    my @results = ();
    my(%seenb4) = ();
    for($j =  $numGenes - 1 ; $j >= 0 ; $j--) { 
	my(@btPath) = ();
        if(defined($seenb4{$bestParses[$j]->[$maxForSt[$j]]})) {
            next;
        }                 
	$seenb4{$maxForSt[$j]} = 1;
	$btPath[$i] = $maxForSt[$i];
	my $len =  abs($sortedPotGenes[$bestParses[$j]->[$maxForSt[$j]]]->{'Exons'}->[0]->[1] - $sortedPotGenes[$bestParses[$j]->[$maxForSt[$j]]]->{'Exons'}->[0]->[2]);
        if($j) {
            for($i = $j ; $i > 0 ; $i--) {
                $btPath[$i - 1] = $btMatrix[$btPath[$i]]->[$i];
		$len += abs($sortedPotGenes[$bestParses[$i - 1]->[$btPath[$i - 1]]]->{'Exons'}->[0]->[1] - $sortedPotGenes[$bestParses[$i - 1]->[$btPath[$i - 1]]]->{'Exons'}->[0]->[2]);
                $seenb4{$bestParses[$j]->[$btPath[$i - 1]]} = 1;
            }
        }
	if(($len * 100.0)/ abs($to - $from) > $minCoverage) {
	    my @genesPath = ();
	    for ($i = 0; $i < scalar(@btPath); $i++) {
		$genesPath[$i] = $sortedPotGenes[$bestParses[$i]->[$btPath[$i]]];
	    }
	    push(@results, [$matScores[$maxForSt[$j]]->[$j], \@genesPath]);
	}
	my @sortedResults = sort Utils::byScore @results;
	return \@sortedResults;
    }
}

sub exonCodingScore {
    my($self, $contig, $lastPos, $nextPos, $strand) = @_;
    my $res = 0.0;
    my $i = 0;
    my $frame = $lastPos % 3;
    if(!defined($strand)) {
	if($lastPos > $nextPos) {
	    $strand = STRAND_COMP;
	    $lastPos = $contig->{'End'} - $lastPos - 1;
	    $nextPos = $contig->{'End'} - $nextPos - 1;
	}
	else {
	    $strand = STRAND_FWD;
	}
    }

    for($i = $lastPos; $i <= $nextPos; $i += 3) {
	$res += $self->{'ContigRegScores'}->{$contig}->{$i}->[$strand]->[EXON];
    }
    return $res*3/($nextPos - $lastPos);
}

sub exonLenScore {
    my($self, $lastTrans, $nextTrans);
    my $len = abs($nextTrans->[TR_POS] - $lastTrans->[TR_POS]);
    return $self->{'HistExonLen'}->smoothedValue($len);
}

sub exonNumScore {
    my($self, $num) = @_;
    return $self->{'HistExonNum'}->smoothedValue($num);
}

sub intronLenScore {
    my($self, $lastTrans, $nextTrans) = @_;
    my $len = abs($nextTrans->[TR_POS] - $lastTrans->[TR_POS]);
    return $self->{'HistIntronLen'}->smoothedValue($len);
}

sub viterbiGene {
    my($self, $bestGenes, $currentExons, $currentScore, $lastTrans, $currNumExons, $numBestGenes, $strand) = @_;
    my $model = $self->{'Model'};
    my($max) = 0;
    my $nextTrans;
    if(defined($lastTrans->[TR_PTR_NEXT]) && scalar(@{$lastTrans->[TR_PTR_NEXT]}) == 0 || !defined($lastTrans->[TR_PTR_NEXT])) {
	return;
    }
#    print join(" | ", @_), "\n";
    foreach $nextTrans (@ { $lastTrans->[TR_PTR_NEXT]}) {
	my $score = $currentScore;
	my $modExon = 0;
	if($nextTrans->[TR_PRED_ST] == EXON) { # comming from a coding state
	    #### only for coding segments
	    my @newExons = ();
	    my $auxScore = 0;
	    my $exon;
	    foreach $exon (@$currentExons) {
		push(@newExons, $exon);
	    }
	    if($lastTrans->[TR_SIG_TYPE] == ACEPTOR) {
		$modExon = 2;
	    }
#	    print ">" x $currNumExons , " ", scalar(@$bestGenes)," $score $nextTrans->[TR_SIG_TYPE] \@ ",$nextTrans->[TR_POS] - 1, " with $lastTrans->[TR_SIG_TYPE] \@ ", $lastTrans->[TR_POS] + $modExon, "\n";
	    
	    push(@newExons, [$nextTrans->[TR_CONTIG], $lastTrans->[TR_POS] + $modExon, $nextTrans->[TR_POS] - 1]);
	    $auxScore = $self->exonCodingScore($lastTrans->[TR_CONTIG], $lastTrans->[TR_POS] + $modExon, $nextTrans->[TR_POS] - 1, $strand);
	    $score += $self->{'Params'}->[CODING]->[$currNumExons]->[$auxScore/$model->{'CodScoreBinLen'}];
	    $score += $self->{'Params'}->[$nextTrans->[TR_SIG_TYPE]]->[$currNumExons]->[($nextTrans->[TR_SIG_SCORE] + 1)/$model->{'SigScoreBinLen'}];
	    $score += $self->{'Params'}->[$lastTrans->[TR_SIG_TYPE]]->[$currNumExons]->[($nextTrans->[TR_SIG_SCORE] + 1)/$model->{'SigScoreBinLen'}];
	    $score += $self->{'Params'}->[EXON_LEN]->[$currNumExons]->[($lastTrans->[TR_POS] - $nextTrans->[TR_POS])/$model->{'HistExonBinLen'}];
	    if($nextTrans->[TR_SIG_TYPE] == STOP) {
		$score += $self->{'Params'}->[EXON_NUM]->[$currNumExons];
		#### We completed a gene, compare it to the rest in the bestGenes array
		if(scalar(@{$bestGenes}) < $numBestGenes) {
		    my $gene = Gene->new('Score' => $score, 'Exons' => \@newExons, 'Strand' => STRAND_FWD);
		    print "New gene ",join(" ", %$gene), " Exons : "; 
		    foreach my $ex (@newExons) { print join(" ",@$ex)," ";}
		    print "\n";
		    $bestGenes->[scalar(@$bestGenes)] = $gene;
		}
		elsif(scalar(@$bestGenes) >=  $numBestGenes && $currentScore > $bestGenes->[$numBestGenes - 1]->{'Score'}) {
		    my $gene = Gene->new('Score' => $score, 'Exons' => \@newExons, 'Strand' => STRAND_FWD);
		    #insert the new gene 
		    print "New gene ",join(" ", %$gene);
		    my $i = scalar(@$bestGenes) - 1;
		    while($i > 0) {
			if($bestGenes->[$i-1]->{'Score'} > $score) {
			    last;
			}
			else {
			    $bestGenes->[$i] = $bestGenes->[$i-1];
			}
			$i--;
		    }
		    $bestGenes->[$i] = $gene;
		}
	    }
	    else {
		##### Recursion #######
		$self->viterbiGene($bestGenes, \@newExons, $score, $nextTrans, $currNumExons + 1, $numBestGenes, $strand);
	    }
	}
	elsif($nextTrans->[TR_PRED_ST] == INTRON) {
	    $score += $self->{'Params'}->[INTRON_LEN]->[$currNumExons]->[($nextTrans->[TR_POS] - $lastTrans->[TR_POS])/$model->{'HistIntronBinLen'}];
	    ##### Recursion #######
	    $self->viterbiGene($bestGenes, $currentExons, $score, $nextTrans, $currNumExons, $numBestGenes, $strand);
	}
	elsif($nextTrans->[TR_PRED_ST] == INTERGENIC) { # a start
	    $self->viterbiGene($bestGenes, $currentExons, $score, $nextTrans, $currNumExons, $numBestGenes, $strand);
	}
	else {
	    die "Dead $nextTrans->[TR_SIG_TYPE]\n";
	}
    }
}


sub addNewTrainingSignals {
    my ($self, $addSignalsFile, $trainingSignals) = @_;
    my $i;
    my $line = "";
    my $seq = "";
    my $header;
    $addSignalsFile ne "" || return;
    open(ADD_SIGNALS, "<$addSignalsFile") || die "$addSignalsFile\n";
    while(defined($header = <ADD_SIGNALS>) && $header =~ /^>(\S+)\s+(\-?\d)/) {
	my $id = $1;
	my $signal = Signal->new('Pos' => $id, 'Class' => $2, 'PredClass' => -1, 'ForTraining' => 1, 'Type' => $i, 'Strand' => STRAND_FWD);
	$seq = <ADD_SIGNALS>;
	chomp($seq);
	$signal->{'Seq'} = $seq;
#	    print join(" ", %$signal),"\n";
	$trainingSignals->{$id} = $signal;
    }
    close(ADD_SIGNALS);
}


sub trainContigSignals {
    my ($self, $org, $addSignals, $minDistToPosSignal, $strand) = @_;
    my $model = $self->{'Model'};
    my $sigType;
    my @svmOpts = (" -t 2 -g 0.1 -c 10 -j 100 ", " -t 2 -g 0.1 -c 10 -j 50 ", " -t 1 -d 2 -c 10 -j 100 ", " -t 2 -g 0.1 -c 10 -j 100 ");
#    for($sigType = START; $sigType <= ACEPTOR; $sigType++) {
    for($sigType = STOP; $sigType <= STOP; $sigType++) {
	print STDERR "Training signal $sigType...\n";
	my %trainingSignals = ();
	my @signals = ();
	$self->filterSignalsForTraining($org, $minDistToPosSignal, $sigType, \%trainingSignals, $strand);
	$self->addNewTrainingSignals($addSignals->[$sigType], \%trainingSignals);
	foreach my $key (keys %trainingSignals) {
	    push(@signals, $trainingSignals{$key});
	}
	&SigUtils::trainSignals(\@signals, $model->{'SignalModels'}->[$sigType], $model->{'SignalFiles'}->[$sigType], $svmOpts[$sigType], $sigType);
    }
}

sub computeIntergenicHistogram {
    my($self, $org, $type) = @_;
    my ($ibeg, $iend);
    my %lenIntergenicVals = ();
    my $class;
    foreach my $contig (keys %{$org->{'Contigs'}}) {
	$ibeg = $org->{'Contigs'}->{$contig}->{'Beg'};
	my @genes = values % { $org->{'Contigs'}->{$contig}->{'Features'}->{$type}};
        $class = $org->{'Contigs'}->{$contig}->{'GCClass'};
        if(!defined($lenIntergenicVals{$class})) {
            %{$lenIntergenicVals{$class}} = ();
        }
	my @sortedGenes = sort {&Gene::compareTo($a, $b);} (@genes);
#	print  scalar(@genes),"\n";
	for(my $i = 0; $i < scalar(@sortedGenes); $i++) {
	    my $gene = $sortedGenes[$i];
	    $iend = $gene->{'Exons'}->[0]->[1] - 1;
	    if($gene->{'Strand'} == STRAND_COMP) {
		$iend = $gene->{'Exons'}->[-1]->[2] - 1;
	    }
	    ($gene->{'Exons'}->[0]->[0] eq $contig) || die;
	    $lenIntergenicVals{$class}->{abs($ibeg - $iend)}++;
	    $ibeg = $gene->{'Exons'}->[-1]->[2] + 1;
	    if($gene->{'Strand'} == STRAND_COMP) {
		$ibeg = $gene->{'Exons'}->[0]->[1] + 1;
	    }	
	}
#	print "mm ", abs($ibeg - $org->{'Contigs'}->{$contig}->{'End'}), "\n";
	$lenIntergenicVals{$class}->{abs($ibeg - $org->{'Contigs'}->{$contig}->{'End'})}++;
    }
    foreach $class (keys %lenIntergenicVals) {
	$self->{'HistIntergenicLen'}->{$class} = Histogram->new('OriginValues' => \%{$lenIntergenicVals{$class}});
	$self->{'HistIntergenicLen'}->{$class}->smooth();
    }
}

sub computeHistograms {
    my($self, $org, $type, $locsByClass) = @_;
    my $numTotal = 0;
    my $lastEnd = -1;
    my $model = $self->{'Model'};
    my $class;
    my $id;
    my %numExonVals = ();
    my %lenExonVals = ();
    my $genes = $org->{'Features'}->{$type};

    foreach $class (keys %$locsByClass) {
	my %lenIntronVals = ();
	foreach $id (@ { $locsByClass->{$class}}) {
	    defined($genes->{$id}) || die "$id does not exist in locations table\n";
	    my $gene = $genes->{$id};
	    my @exons = @ { $gene->{'Exons'}};
	    my $i;
	    $lastEnd = 0;
	    $numExonVals{scalar(@exons)}++;
	    $numTotal++;
	    for($i = 0; $i < scalar(@exons);$i++) {
		$lenExonVals{abs($exons[$i]->[2] - $exons[$i]->[1]) + 1}++;
		if($lastEnd > 0) {
		    $lenIntronVals{abs($exons[$i]->[1] - $lastEnd)}++;
		}
		$lastEnd = $exons[$i]->[2];
	    }
	}
	$self->{'HistIntronLen'}->{$class} = Histogram->new('OriginValues' => \%lenIntronVals);
	$self->{'HistIntronLen'}->{$class}->smooth();
#	print $class, " ", scalar(keys %lenIntronVals), "\n";
#	print $class, "INTRONLEN ", scalar(keys %{$self->{'HistIntronLen'}->{$class}->{'SmoothValues'}}), "\n";

    }

    $self->{'HistExonNum'} = Histogram->new('OriginValues' => \%numExonVals); 
    $self->{'HistExonNum'}->smooth();
    $self->{'HistExonLen'} = Histogram->new('OriginValues' => \%lenExonVals);
    $self->{'HistExonLen'}->smooth();
#    print scalar(keys %numExonVals), " ", scalar(keys %lenExonVals) ,"\n";
#    print "EXONLEN ", scalar(keys %{$self->{'HistExonNum'}->{'SmoothValues'}}), " ", scalar(keys %{$self->{'HistExonLen'}->{'SmoothValues'}}), "\n";
# Go for intergenic regions
    $self->computeIntergenicHistogram($org, $type);
    ### returning the three histograms obtained from this set of sequences ###
    return [$self->{'HistExonLen'}, $self->{'HistExonNum'}, $self->{'HistIntronLen'}, $self->{'HistIntergenicLen'}];
}
    
sub trainCodingRegions {
    my($self, $genes, $locsByClass) = @_;
    my %empty = ();
    my $model = $self->{'Model'};
    my $class;
    my $file = $model->{'CodingModel'}.".fsa";    
    open(FILE, ">$file") || die "Can't open $file\n";     
    foreach $class (keys %$locsByClass) { 
	my $id; 
	foreach $id (@ {$locsByClass->{$class}}) { 
	    my $gene = $genes->{$id};
 	    print FILE ">",$gene->{'Id'},"\n", $gene->{'Seq'},"\n"; 
	}
    }
    close(FILE);
    system("gfinder_train -d $model->{'Domain'} -printGC $model->{'GCProp'} -t $file -R -M $model->{'CodingModel'} > /dev/null"); 
#	unlink($file); 
    foreach $class (keys %$locsByClass) { 
	`cp $model->{'CodingModel'} "$model->{'CodingModel'}.$class"`;
    } 
}

sub trainIStateRegions {
    my($self, $org, $type) = @_;
    my %empty = ();
    my $model = $self->{'Model'};
    my $class;

    my $fileIntron = $model->{'IntronModel'}.".fsa";
    my $fileIntergenic = $model->{'IntergenicModel'}.".fsa";
    open(FILE, ">$fileIntron") || die "Can't open $fileIntron\n"; 
    open(FILEIG, ">$fileIntergenic") || die "Can't open $fileIntergenic\n"; 

    my $contigs = $org->{'Contigs'};

    for my $c (keys %$contigs) {
        my @array = ();
        my $i;
        my $gene;
        my $strand;
        for($i = 0; $i <= $contigs->{$c}->{'Len'}; $i++) {
            push(@array, 0);
        }
        my @genes = sort {&Gene::compareTo($a, $b);} (values %{ $contigs->{$c}->{'Features'}->{'TRAIN'}});
#    print "$c ", $contigs->{$c}->{'Len'}," ", scalar(@genes), "\n";
        
        my $lastGeneId = "";
        my $lastMin = 1;
        my $geneId = "";

        for $gene (@genes) {
            my $min = 0;
            my $max = 0;
            $strand = $gene->{'Strand'};
            if($gene->{'Id'} =~ /(\S+)\|([^\-]+)$/) {
                $geneId = $1;
            }
            else {
                $geneId = $gene->{'Id'};
            }
            
            if($lastGeneId ne "" && $geneId eq $lastGeneId) {
                $min = $lastMin;
            }
            else {
                $min = $gene->lowestCoord();
            }
            
            $max = $gene->highestCoord();
            
#        print STDERR "$gene->{'Id'} \"$min\" \"$max\" $contigs->{$c}->{'Len'}\n";
            for($i = $min; $i <= $max; $i++) {
                $array[$i] = 1;
            }

            $lastGeneId = $geneId;
            $lastMin = $min;;
            
            my $exons = $gene->{'Exons'};
            
            if(scalar(@{$exons}) == 1) {
                next;
            }

            for(my $i = 1; $i < scalar(@{$exons}); $i++) {
                my $beg = $exons->[$i - 1]->[2] + 1;
                my $end = $exons->[$i]->[1] - 1;
                if(abs($end - $beg) < 6) {
                    next;
                }
                if($gene->{'Strand'} == STRAND_COMP) {
                    $beg -= 2;
                    $end += 2;
                }
                print FILE ">",$gene->{'Id'},"_".$beg."_"."$end\n";
                print FILE $contigs->{$c}->getSequence(0,$beg, $end,0, 0, $gene->{'Strand'}),"\n"; 
            }
            
        }
        
        my $interBeg = -1;
        my $interEnd = -1;
        for($i = 1; $i <= $contigs->{$c}->{'Len'}; $i++) {
            if($array[$i] == 0 && $interBeg == -1) {
                $interBeg = $i;
            }
            if($i > 1) {
                if(($array[$i] == 1  || $i == $contigs->{$c}->{'Len'}) && $array[$i - 1] == 0) {
                    $interEnd = $i;
                }
            }

            if($interBeg != -1 && $interEnd != -1) {
                if($interEnd - $interBeg  < 6) {
                    $interBeg = -1;
                    $interEnd = -1;
                    next;
                }
                if($strand == STRAND_COMP) {
                    my $tmp = $interBeg;
                    $interBeg = $interEnd;
                    $interEnd = $tmp;
                }
                
                print FILEIG ">",$c,"_".$interBeg."_"."$interEnd\n";
                print FILEIG $contigs->{$c}->getSequence(0,$interBeg, $interEnd, 0, 0, $strand),"\n"; 
                $interBeg = -1;
                $interEnd = -1;
                
            }        
        }
    }
    
    close(FILE);
    close(FILEIG);
#    system("istate_train -d $model->{'Domain'} -t $fileIntron -R -M $model->{'IntronModel'} > /dev/null"); 
#    system("istate_train -d $model->{'Domain'} -t $fileIntergenic -R -M $model->{'IntergenicModel'} > /dev/null"); 
#   unlink($fileIntron); 
#   unlink($fileIntergenic); 
}

sub genGCIsochores {
    my($self, $contigs, $gcClasses, $scoreDir) = @_;
    open(CL, ">$scoreDir/clusters") || die "Can't open $scoreDir/clusters for writing\n";
    foreach my $c (keys %$contigs) {
	my $contig = $contigs->{$c};
	my $gcClass = SigUtils::getGCClass($contig->{'Seq'}, $gcClasses);
	$contig->{'GCClass'} = $gcClass;
	print CL "$contig->{'Id'} ", STRAND_FWD, " $contig->{'Beg'} ", length($contig->{'Seq'}), " $gcClass\n";
	print CL "$contig->{'Id'} ", STRAND_COMP, " ", length($contig->{'Seq'}), " $contig->{'Beg'} $gcClass\n";
    }
    close(CL);
}

# this routine is useful when contigs are extremely long, but only one isochore
# which pretty much covers 99% of cases. Used mostly for testing
sub scoreGCIsochores {
    my($self, $contigs, $gcClasses, $scoreDir) = @_;
    foreach my $c (keys %$contigs) {
	my $contig = $contigs->{$c};
	$contig->scoreContig($self->{'Model'}, $gcClasses, $scoreDir); 
    }
}


sub genContigClusters {
    my($self, $contigs, $gcClasses, $scoreFile) = @_;    
    # getting the length of the possible gene in the training set
    my $window = 1000000;
    my $contig;
    my $model = $self->{'Model'};
    my %contigClusters = %{$self->{'ContigGCClusters'}} = %{&SigUtils::getContigClustersBySignal($contigs, $gcClasses, $window, "G|C", 1)};
    my $tmpFile = "$model->{'CodingModel'}.coding";

    foreach my $gcClass (%$gcClasses) {
	if(-e "$scoreFile.in.$gcClass") {
	    unlink "$scoreFile.in.$gcClass";
	}
    }
    # preparing input files
    foreach $contig (keys %contigClusters) {
	my $strand; 
	for($strand = STRAND_FWD ; $strand <= STRAND_COMP ; $strand++) {
	    my $cluster;
	    # flushing out all the contig regions that belong to this class
	    foreach $cluster (@{$contigClusters{$contig}->[$strand]}) {
		open(IN, ">>$scoreFile.in.$cluster->[GC_CLUSTER_CLASS]") || die "Can't open $scoreFile.in.$cluster->[GC_CLUSTER_CLASS]\n";
		print IN ">$contig\:$strand\:", $cluster->[GC_CLUSTER_BEG] - $contigs->{$contig}->{'Beg'} + 1,"\n";
		print IN $contigs->{$contig}->getSequence(0,$cluster->[GC_CLUSTER_BEG],$cluster->[GC_CLUSTER_END], 0, 1, $strand),"\n";
		close(IN);
	    }
	}
    }
}

sub scoreContigClusters {
    my($self, $gcClasses, $scoreFile) = @_;    
    # getting the length of the possible gene in the training set
    my $model = $self->{'Model'};
    my $command = "";
    foreach my $gcClass (keys %$gcClasses) {
	my $thisModel = "$model->{'CodingModel'}.$gcClass"; 

	(-e "$scoreFile.in.$gcClass") || next;
	if(! -e "$model->{'CodingModel'}.$gcClass" || -z "$model->{'CodingModel'}.$gcClass") {
	    print STDERR "Warning : there weren't enough training sequences for the $gcClass GC class. Model hasn't been computed, using alternative model\n";
	    my $found = 0;
	    foreach my $gcClassAlt (keys %$gcClasses) {
		if(-e "$model->{'CodingModel'}.$gcClassAlt" && !-z $model->{'CodingModel'}.".$gcClassAlt") {
		    $thisModel = "$model->{'CodingModel'}.$gcClassAlt";
		    $found = 1;
		    last;
		}
	    }
	    ($found == 1) || die "Can't find coding model\n";
	}
	$command = "score_contigs -d $model->{'Domain'} -f $scoreFile.in.$gcClass -M $thisModel > $scoreFile.immcoding.$gcClass";
	system($command);    
#	$command = "score_contigs_nn -d $model->{'Domain'} -f $scoreFile.in.$gcClass > $scoreFile.coding.$gcClass";
#	system($command);    
	$command = "score_contigs -d $model->{'Domain'} -f $scoreFile.in.$gcClass -M $model->{'IntronModel'} > $scoreFile.intron.$gcClass";
	system($command);
	$command = "score_contigs -d $model->{'Domain'} -f $scoreFile.in.$gcClass -M $model->{'IntergenicModel'} > $scoreFile.intergenic.$gcClass";
	system($command);
    }
}

sub storeClusters {
    my($self, $file, $contigs) = @_;
    my %contigClusters = %{$self->{'ContigGCClusters'}};
    my $contig;
    
    open(CL, ">$file.clusters") || die "Can't open $file.clusters for writing\n";
    foreach $contig (keys %contigClusters) {
	my $strand;
	for($strand = STRAND_FWD ; $strand <= STRAND_COMP ; $strand++) {
	    foreach my $cluster (@{$contigClusters{$contig}->[$strand]}) {
		print CL "$contig $strand ", $cluster->[GC_CLUSTER_BEG] - $contigs->{$contig}->{'Beg'} + 1, " ", $cluster->[GC_CLUSTER_END] - $contigs->{$contig}->{'Beg'} + 1, " $cluster->[GC_CLUSTER_CLASS]\n";
	    }
	}
    }
    close(CL);
}


sub retrieveContigRegScoresFromFile {
    my($self, $file, $contigs, $type) = @_;
    my $line = <$file>;
    my $i;
    my ($contig, $strand, $begin) = (0,0,0);
    while(defined($line) && $line =~ /^\>(\S+)\:(\d)\:(\d+)/) {
	($contig, $strand, $begin) = ($1,$2,$3);
	$i = 0;
	while(defined($line = <$file>) && $line !~ /^\>/) {
	    chomp($line);
	    $begin = $begin + $contigs->{$contig}->{'Beg'} - 1;
	    if($strand == STRAND_FWD)  {
		$self->{'ContigRegScores'}->{$contig}->{$begin + $i++}->[$strand]->[$type] = $line;
#		print STDERR "$type $contig $begin $i $strand $self->{'ContigRegScores'}->{$contig}->{$begin + $i - 1}->[$strand]->[$type]\n";
	    }
	    else {
		
		$self->{'ContigRegScores'}->{$contig}->{$begin - $i++}->[$strand]->[$type] = $line;
	    }
	}
    }
}


sub retrieveContigRegScores {
    my($self, $scoreFile, $gcClasses, $contigs) = @_;
    my $line = "";
    foreach my $gcClass (keys %$gcClasses) {
#	print STDERR "$scoreFile.coding.$gcClass", "\n";
	open(EXN, "<$scoreFile.immcoding.$gcClass") || next;
	open(INTR, "<$scoreFile.intron.$gcClass") || next;
	open(INTE, "<$scoreFile.intergenic.$gcClass") || next;
	$self->retrieveContigRegScoresFromFile(\*EXN, $contigs, EXON);
	$self->retrieveContigRegScoresFromFile(\*INTR, $contigs, INTRON);	
	$self->retrieveContigRegScoresFromFile(\*INTE, $contigs, INTERGENIC);	
	close(EXN);
	close(INTR);
	close(INTE);
    }
    
    open(CL, "<$scoreFile.clusters") || die "Can't open $scoreFile.clusters for reading\n";
    while(defined($line = <CL>)) {
	chomp($line);
	my @line = split(" ", $line);
	$line[2] = $line[2] + $contigs->{$line[0]}->{'Beg'} - 1;
	$line[3] = $line[3] + $contigs->{$line[0]}->{'Beg'} - 1;
	push(@{$self->{'ContigGCClusters'}->{$line[0]}->[$line[1]]}, [$line[2], $line[3], $line[4]]);
    }
    close(CL);
}


sub scoreContigReg {
    my ($self, $contig, $beg, $end, $type, $frame) = @_;
    my $i;
    my $score = 0;
    my $step = 3;
    $frame = 0 unless defined($frame);
    if($beg + $frame < 1) {
	$beg = 1 - $frame;
    }
    if($end < 1) {
	$end = 1;
    }
    if($beg + $frame > $contig->{'Len'}) {
	$beg = $contig->{'Len'} - $frame;
    }
    if($end > $contig->{'Len'}) {
	$end = $contig->{'Len'};	
    }
    if($type != EXON) {
	$step = 1;
    }
    if($end < $beg) {
	for($i = $beg + $frame; $i > $end; $i -= $step) {
	    my $sc = $self->{'ContigRegScores'}->{$contig->{'Id'}}->{$i}->[STRAND_COMP]->[$type];
	    for(my $j = 0; $j < $step; $j++) {
		$score += log($sc + 1e-50)/log(2.0);
	    }
	}
    }
    else {
#	print STDERR "$contig->{'Id'} $beg $end $type\n";
	for($i = $beg + $frame; $i < $end; $i += $step) {
#	    print STDERR $i, " ",  $self->{'ContigRegScores'}->{$contig->{'Id'}}->{$i}->[STRAND_FWD]->[$type], " ", $score, "\n";
	    my $sc = $self->{'ContigRegScores'}->{$contig->{'Id'}}->{$i}->[STRAND_FWD]->[$type];
	    for(my $j = 0; $j < $step; $j++) {
#		$score += log($sc + 1e-50)/log(2.0);
		$score += $sc;
	    }
	}
    }
    return $score;
}


sub scoreContigSignals {
    my($self, $org, $maxNumExons, $verbose, $scoreFile, $tag) = @_;
    my $model = $self->{'Model'};
    my $i = 0;
    my $j = 0;
    my $index;
    my $accIndex = 0;
    $tag = "" if !defined($tag);
    for($i = START; $i <= ACEPTOR;$i++) {
#    for($i = DONOR ; $i <= DONOR ; $i++) {
	$accIndex = 0;
	$index = 0;
	my @signals = undef;
	my @testSignals = ();
	my $signal = undef;
	if($verbose) {
	    print STDERR "Finding signal \#", $i, " in contigs\n";
	}
	foreach my $cKey (keys % {$self->{'ContigSigScores'}}) {
	    my $contig = $org->{'Contigs'}->{$cKey};
	    foreach my $ind (keys %{$self->{'ContigSigScores'}->{$cKey}}) {
		for(my $strand = STRAND_FWD; $strand <= STRAND_COMP;$strand++) {
		    if(!defined($self->{'ContigSigScores'}->{$cKey}->{$ind}->[$strand]) || ($self->{'ContigSigScores'}->{$cKey}->{$ind}->[$strand]->[0]&~COMP) != $i) {
			next;
		    }
		    $signal = Signal->new('Pos' => $cKey."_".$ind, 'Class' => $self->{'ContigSigScores'}->{$cKey}->{$ind}->[$strand]->[2], 'PredClass' => -1, 'ForTraining' => 1, 'Type' => $self->{'ContigSigScores'}->{$cKey}->{$ind}->[$strand]->[0], 'Strand' => $strand);
		    my $numExons = $maxNumExons;
		    my $exon = $signal->getAssocExon($numExons, int(rand $numExons));
		    #print join(" ",%$signal),"\n";
		    $signal->{'ForTraining'} = 1;
		    $signal->initSeq($org->{'Contigs'}, $model->{'SignalWindows'}->[$i]);
		    if($signal->{'ForTraining'} == 0) {
			$self->{'ContigSigScores'}->{$cKey}->{$ind}->[$strand]->[1] = -2.0;
			next;
		    }
		    push(@testSignals, $signal);
		    $signals[$index] = $self->{'ContigSigScores'}->{$cKey}->{$ind}->[$strand];
		    if(!($index++ % 10000)) {
			print STDERR ".";
		    }
		    $accIndex++;
		    if(!($index % 750000)) {
			my @results = @ { &SigUtils::testSignals(\@testSignals, $model->{'SignalModels'}->[$i], $self->{'SignalFiles'}->[$i], $i)};
			for(my $r = 0 ; $r < scalar(@results); $r++) {
			    $signals[$r]->[1] = $results[$r];
			}
			@testSignals = ();
			@signals = undef;
			$index = 0;
		    }
		}
	    }
	}
	if($verbose) {
	    print STDERR " Total number so far is ", $accIndex, "... Scoring them .. \n";
	}
	if(!$index) { # no more signals were added
	    next;
	}
	my @results = @ { &SigUtils::testSignals(\@testSignals, $model->{'SignalModels'}->[$i], $model->{'SignalFiles'}->[$i], $i, $tag)};
	for(my $r = 0 ; $r < scalar(@results); $r++) {
	    $signals[$r]->[1] = $results[$r];
	}
    }
}   

sub storeContigSigScores { 
    my($self, $file, $contigs, $tag) = @_;
    $tag = "" if !defined($tag);
    open(TMP, ">$file"."$tag.signals") || die "Can't open $file"."$tag.signals for writing\n";

    foreach my $c (sort {$a cmp $b} keys %{$self->{'ContigSigScores'}}) {
	foreach my $p (sort {$a <=> $b} keys %{$self->{'ContigSigScores'}->{$c}}) {
#	    print $c, " ", $p," ",$p - $contigs->{$c}->{'Beg'} + 1," ";
#	    print $self->{'ContigSigScores'}->{$c}->{$p}->[STRAND_COMP]->[0], " - ";
#	    print $self->{'ContigSigScores'}->{$c}->{$p}->[STRAND_FWD]->[0], " - ";
	    for(my $strand = STRAND_COMP; $strand >= STRAND_FWD; $strand--) {
#		print $self->{'ContigSigScores'}->{$c}->{$p}->[$strand]->[0]," + ";
		my $position = $p;
		if(!defined($self->{'ContigSigScores'}->{$c}->{$p}->[$strand]->[0])){
#		    print "bad\n";
		    next;
		}
#		print "good\n";

		my @content = @{$self->{'ContigSigScores'}->{$c}->{$p}->[$strand]};
		my $oldPos = $position;
		$position = $position - $contigs->{$c}->{'Beg'} + 1;
		if($strand == STRAND_COMP) {
		    $position = $contigs->{$c}->{'Len'} - $position + 1;
		}
		print TMP "$c ", $oldPos, " ", $position, " $strand ", join(" ",@content),"\n";
	    }	    
	}
    }
    close(TMP);
}

sub retrieveContigSigScores {
    my($self, $file, $contigs, $sigType) = @_;
    open(TMP, "<$file.signals") || die "Can't open $file.signals for reading\n";
    my $line = "";
    my @line = ();
    while(defined($line = <TMP>)) {
	chomp($line);
	@line = split(" ", $line);
	if(defined($sigType) && $sigType != $line[3]) {
	    next;
	} 
	if($line[3] == STRAND_COMP) {
	    $line[2] = $contigs->{$line[0]}->{'Len'} - $line[2] + 1;
	}
	$line[2] = $line[2] + $contigs->{$line[0]}->{'Beg'} - 1;
	($line[1] == $line[2]) || die "Signal entry at ".join(" ", @line)."\n";
	$self->{'ContigSigScores'}->{$line[0]}->{$line[2]}->[$line[3]] = [$line[4], $line[5], $line[6]];
    }
    close(TMP);
}

sub setContigSignals {
    my ($self, $signals) = @_;
    $self->{'ContigSigScores'} = $signals;
}

sub classifySignalsWGeneSet {
    my($self, $org, $type, $minSepSignals) = @_;
    my ($id, $i, $shift) = (0,0,0);
    defined(my $genes = $org->{'Features'}->{$type}) || die "Orf set $type not specified\n";
    my $contigs = $org->{'Contigs'};
    my $signals = $self->{'ContigSigScores'};
    my $signal;
    my $window = undef;
    foreach $id (keys %$genes) {
	my $gene = $genes->{$id};
        # checking contig id in the first exon
        if(!defined($signals->{$gene->{'Exons'}->[0]->[0]})) {
            next;
        }
        defined($contigs->{$gene->{'Exons'}->[0]->[0]}) || die "contig $gene->{'Exons'}->[$i]->[0]}\n";
	for($i = 0; $i < scalar(@ {$gene->{'Exons'}}); $i++) {
            my $exon = $gene->{'Exons'}->[$i];
            my $signalInd;
            if($i == 0) { # getting the start
		if($gene->{'Strand'} == STRAND_FWD) {
                    $shift = 0;
                    $signalInd = START;
		}
                else { # complementary signal
                    $shift = 0;
                    $signalInd = START|COMP;
		}
            }
            else { # getting aceptor
		if($gene->{'Strand'} == STRAND_FWD) {
                    $signalInd = ACEPTOR;
                    $shift = -2;
		}
                else {
                    $signalInd = ACEPTOR|COMP;
                    $shift = 2;
		}
            }

	    $signal = $signals->{$exon->[0]}->{$exon->[1] + $shift}->[$gene->{'Strand'}];
#	    print "Calling filter with ", join(" ", @$exon),"\n";

	    if(!$gene->{'NoStart'} || ($signalInd&~COMP) != START) { 
		if(defined($signal) && ($signal->[0] == ($signalInd&~COMP))) {
		    $signal->[2] = 1;
		}
		else {
                    if($exon->[1] + $shift >= 1 && $exon->[1] + $shift <= $contigs->{$gene->{'Exons'}->[0]->[0]}->{'Len'}) {
                        print STDERR "Warning : signal ".$signalInd." is not specified at $exon->[0] $exon->[1] + $shift ...\n";
                        $signals->{$exon->[0]}->{$exon->[1] + $shift}->[$gene->{'Strand'}] = [$signalInd&~COMP, 0, 1];
                    }
		}
	    }

            if($i == scalar(@ {$gene->{'Exons'}}) - 1) { # getting the stop
		if($gene->{'Strand'} == STRAND_FWD) {
                    $shift = 1;
                    $signalInd = STOP;
		}
                else {
                    $shift = -1;
                    $signalInd = STOP|COMP;
		}
            }
            else { # getting  donors
		if($gene->{'Strand'} == STRAND_FWD) {
                    $shift = 1;
                    $signalInd = DONOR;
		}
                else {
                    $shift = -1;
                    $signalInd = DONOR|COMP;
		}
            }
	    $signal  = $signals->{$exon->[0]}->{$exon->[2] + $shift}->[$gene->{'Strand'}];
#	    print "Calling filter with ", join(" ", @$exon),"\n";

	    if(!$gene->{'NoStop'} || ($signalInd&~COMP) != STOP) { 
		if(defined($signal) && ($signal->[0] == ($signalInd&~COMP))) {
		    $signal->[2] = 1;	    
		}
		else {
                    if($exon->[2] + $shift >= 1 && $exon->[2] + $shift <= $contigs->{$gene->{'Exons'}->[0]->[0]}->{'Len'}) {
                        $signals->{$exon->[0]}->{$exon->[2] + $shift}->[$gene->{'Strand'}] = [$signalInd&~COMP, 0, 1];
                        print STDERR "Warning : signal ".$signalInd." is not specified at $exon->[0] $exon->[2] + $shift...\n";
                    }
		}
	    }
	}
    }
}


sub filterSignalsForTraining {
    my($self, $org, $minDistToPosSig, $sigType, $trainingSignals, $str, $noSeq) = @_;
    my ($id, $i, $shift) = ("",0,0);
    my $contigs = $org->{'Contigs'};
    my $signals = $self->{'ContigSigScores'};
    my $signal;
    my @strands = ();
    $noSeq = 0 if !defined($noSeq);
    if(!defined($str)) {
	@strands = (STRAND_FWD, STRAND_COMP);
    }
    else {
	@strands = ($str);
    }
    # #initializing contig area
    foreach my $cKey (keys %$contigs) {
	my @posSignals = sort {$a <=> $b} keys %{$signals->{$cKey}};
	my @lastInd = (-1,-1);
	my @lastPos = (0,0);
	for($i = 0; $i < scalar(@posSignals); $i++) {
	    foreach my $strand (@strands) {
                my $position = $posSignals[$i];
		my $sig = $signals->{$cKey}->{$position}->[$strand];
                if(!defined($sig)) {
#                   print "bad\n";
		    next;
		}
		if($sig->[0] != $sigType) {
		    next;
		}

		if($lastInd[$strand] == -1 || $position >= $lastPos[$strand] + $minDistToPosSig) {
		    my $p = $lastInd[$strand] + 1;
		    while($p < scalar(@posSignals)) {
			if(defined($signals->{$cKey}->{$posSignals[$p]}->[$strand])) {
			    if($sig->[0] == $signals->{$cKey}->{$posSignals[$p]}->[$strand]->[0] && $signals->{$cKey}->{$posSignals[$p]}->[$strand]->[2] == 1) {
				last;
			    }
			}
			$p++;
		    }
		    $lastInd[$strand] = $p;
		    if($lastInd[$strand] >= scalar(@posSignals)) {
			$lastPos[$strand] = $contigs->{$cKey}->{'Len'} + $minDistToPosSig;
		    }
		    else {
			$lastPos[$strand] = $posSignals[$lastInd[$strand]];
		    }
		}

		if($sig->[2] == 1 || $position < $lastPos[$strand] - $minDistToPosSig) {
		    my $sigLoc = $cKey."_".$position;
		    $signal = Signal->new('Pos' => $sigLoc, 'Class' => $sig->[2], 'PredClass' => -1, 'ForTraining' => 1, 'Type' => $sig->[0], 'Strand' => $strand, 'Score' => $sig->[1]);
		    if(!$noSeq) {
			$signal->initSeq($contigs, $self->{'Model'}->{'SignalWindows'}->[$sig->[0]&~COMP]);
		    }
#		    print $lastPos[$strand]," ", $lastInd[$strand], " " ,join(" ", %$signal), "\n";

		    if($signal->{'ForTraining'} == 1) {
#	    print "Added signal ", join(" ", %$signal),"\n";
			$trainingSignals->{$sigLoc} = $signal;
		    }
		    elsif($signal->{'Class'} == 1) {
			print "Bad Signal at $sigLoc : $sig->[0]\n";  
		    }
		}
	    }
	}
    }
}    


# Filter signals from the set of all possible ones, given a positive signal. MinSepSignals and ratioPosNeg are used for this. Negative examples will be chosen from the ones around the positive one
sub filterAroundPosSignal {
    my($self, $contigSignals, $window, $pSContig, $pSCoord, $pSStrand, $upsLimit, $minDistToPosSig, $contigs, $trainingSignals, $upsOnly) = @_;
    
    my $sigLoc = $pSContig."_".$pSCoord;
#    print join(" ", @{$gene->{'Exons'}->[$i]}), " $sigLoc ", $i, " ", $signalInd,"\n";
    my $sig = $contigSignals->{$pSCoord}->[$pSStrand];
    my $signalInd = $sig->[0];
    my $signal = Signal->new('Pos' => $sigLoc, 'Class' => 1, 'PredClass' => -1, 'ForTraining' => 1, 'Type' => $sig->[0], 'Strand' => $pSStrand);
    my $signalUps = $signal;
    my $stepUps = -1;
    my $lastUps = $pSCoord  + $stepUps*$minDistToPosSig;
    my $stepDws = 1;
    my $lastDws = $pSCoord  + $stepDws*$minDistToPosSig;
    
    if(!defined($trainingSignals->[$sig->[0]])) {
	%{$trainingSignals->[$sig->[0]]} = ();
    }
    
    #generating the negative examples
    my $insert = 1; 
    while(defined($signal)) {
	#updating signal information
#	print "Init seq $window\n";

	$signal->{'ForTraining'} = 1;
	$signal->initSeq($contigs, $window);
	if($signal->{'ForTraining'} == 1 && $insert == 1) {
#	    print "Added signal ", join(" ", %$signal),"\n";
	    $trainingSignals->[$sig->[0]]->{$sigLoc} = $signal;
	}
	else {
	    if($insert == 1 && $signal->{'Class'} == 1) {
		print "Bad Signal at $sigLoc : $sig->[0]\n";
	    }
	}
	# alternating signals from upstream and downstream to find negative signals
	$insert = 0;
	if(defined($signalUps)) {
            do {
                $lastUps += $stepUps;
		$sigLoc = $pSContig."_".$lastUps;
#		print "Ups.. $sigLoc $gene->{'Strand'} ", defined($contigSignals->{$lastUps}->[$gene->{'Strand'}]), "\n";
            } while(!defined($contigSignals->{$lastUps}->[$pSStrand]) && $lastUps >= $contigs->{$pSContig}->{'Beg'});
#	    print defined($contigSignals->{$lastUps}->[$gene->{'Strand'}]), "\n";
            if(defined($contigSignals->{$lastUps}->[$pSStrand]) && !defined($trainingSignals->[$signalInd]->{$sigLoc}) && $lastUps >= $upsLimit) {
#		print "Getting signal type $contigSignals->{$lastUps}->[$gene->{'Strand'}]->[0]\n";
		if(!defined($trainingSignals->[$signalInd]->{$sigLoc}) && $contigSignals->{$lastUps}->[$pSStrand]->[0] == $signalInd) {
		    $signalUps = $signal = Signal->new('Pos' => $sigLoc, 'Class' => -1, 'PredClass' => -1, 'ForTraining' => 1, 'Type' => $signalInd, 'Strand' => $signalInd&COMP);
#		    print join(" ", %$signal),"\n";
		    $insert = 1;
		}
            }
	    else {
                $signalUps = $signal = undef;
            }
	}
        elsif(!$upsOnly) { # run out of signals going ups... trying going ds
            do {
                $lastDws += $stepDws;
		$sigLoc = $pSContig."_".$lastDws;
#		print "Down.. ", $sigLoc," $gene->{'Strand'} ", defined($contigSignals->{$lastDws}->[$gene->{'Strand'}]),"\n";
            }while(!defined($contigSignals->{$lastDws}->[$pSStrand]) && $lastDws < $contigs->{$pSContig}->{'Beg'} + $contigs->{$pSContig}->{'Len'});
#	    print defined($contigSignals->{$lastDws}->[$gene->{'Strand'}]), "\n";
	    if(defined($contigSignals->{$lastDws}->[$pSStrand]) && $lastDws <= $contigs->{$pSContig}->{'End'}) {
		if(!defined($trainingSignals->[$lastDws]->{$sigLoc}) && $contigSignals->{$lastDws}->[$pSStrand]->[0] == $signalInd) {
#		    print "Getting signal type $contigSignals->{$lastDws}->[$gene->{'Strand'}]->[0]\n";
		    $signal = Signal->new('Pos' => $sigLoc, 'Class' => -1, 'PredClass' => -1, 'ForTraining' => 1, 'Type' => $signalInd, 'Strand' => $signalInd&COMP);
#		    print join(" ", %$signal),"\n";
		    $insert = 1;
		}
	    }
	    else {
		$signal = undef;
	    }
	}
    }
}
1;
