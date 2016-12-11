package Gene;
require Exporter;
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Feature;
no warnings 'redefine';
@ISA = qw (Feature)
;@EXPORT= qw (
              moveLocation2OrigContigs
              getFeatureSeq
              shiftCoords
              lowestCoord
              highestCoord
              equalTo
              compareTo
              displayHeader
              getThreePWeight
              getFivePWeight
              getCDSWeight
              dump
              dumpUTRSeqs
              cropCoords
              contigId
              addUTRAnnot
              );
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
    
sub new {
    my $type = shift;
    my %params = @_;

    if(defined($params{'Gene'})) {
        my %gene = %{$params{'Gene'}};
        
        my $i;
        my $j;
        %params = %gene;

        undef $params{'Exons'};
        @{$params{'Exons'}} = ();
        for($i = 0; $i < scalar(@{$gene{'Exons'}}); $i++) {
            @{$params{'Exons'}->[$i]} = @{$gene{'Exons'}->[$i]};
        }
        undef $params{'5UTRExons'};
        @{$params{'5UTRExons'}} = ();
        for($i = 0; $i < scalar(@{$gene{'5UTRExons'}}); $i++) {
            @{$params{'5UTRExons'}->[$i]} = @{$gene{'5UTRExons'}->[$i]};
        }
        undef $params{'3UTRExons'};
        @{$params{'3UTRExons'}} = ();
        for($i = 0; $i < scalar(@{$gene{'3UTRExons'}}); $i++) {
            @{$params{'3UTRExons'}->[$i]} = @{$gene{'3UTRExons'}->[$i]};
        }

    }
    
    @{$params{'Exons'}} = () if !defined($params{'Exons'});
    
    my $self = Feature::new($type, %params);

    $params{'Phase'} = 0 if !defined($params{'Phase'});
    $params{'NoStart'} = 0 if !defined($params{'NoStart'});
    $params{'NoStop'} = 0 if !defined($params{'NoStop'});
    @{$params{'5UTRExons'}} = () if !defined($params{'5UTRExons'});
    @{$params{'3UTRExons'}} = () if !defined($params{'3UTRExons'});

    $self->{'Phase'} = $params{'Phase'};
    $self->{'NoStart'} = $params{'NoStart'};
    $self->{'NoStop'} = $params{'NoStop'};
    $self->{'gId'} = $params{'gId'};
    $self->{'tId'} = $params{'tId'};
    
    if(!defined($params{'Id'})) {
        $params{'Id'} = $params{'gId'};
        if(defined($params{'tId'})) {
            $params{'Id'} .= "|".$params{'tId'};
        }
    }
    
    $self->{'Id'} = $params{'Id'};
    
    $self->{'5UTRExons'} = $params{'5UTRExons'};
    $self->{'3UTRExons'} = $params{'3UTRExons'};

    ($self->{'NoStart'} != 0 || ($self->{'Phase'} == 0 && $self->{'NoStart'} == 0)) || warn "Phase can't be ".$self->{'Phase'}." if start is present for gene $self->{'Id'}\n";
    $self->getFeatureSeq();
#    print "contig ", $self->{'Contig'}," ", $self->{'Seq'}, "\n";
    return $self;
    
}

sub checkStart {
    my($self) = @_;
    my $frame =  (3 - $self->{'Phase'}) % 3;
    if($frame) {
	$self->{'NoStart'} = 1;
	return;
    }

    if($self->{'NoStart'} && substr($self->{'Seq'}, 0, 3) =~ /atg/i) {
        my $numStops = GenUtils::numPhaseOccurrences($self->{'Seq'}, 0, "taa", "tga", "tag");
        if(!$numStops) {
	    $self->{'Phase'} = 0;
	    $self->{'NoStart'} = 0;
	}
    }

    if(!$self->{'NoStart'} && substr($self->{'Seq'}, 0, 3) !~ /atg/i) {
        $self->{'NoStart'} = 1;
    }
}

sub checkStop {
    my($self) = @_;
    my $frame =  (3 - $self->{'Phase'}) % 3;
    my $endPhase = (length($self->{'Seq'}) - $frame) % 3;
#    print STDERR "phases ", $frame, " ", length($self->{'Seq'}), " ", $self->{'Phase'}, " ", $endPhase, "\n";
    my $stop1 = $endPhase ? "nnn" : substr($self->{'Seq'}, length($self->{'Seq'}) - 3, 3);
    my $stop2 = $endPhase ? "nnn" : $self->{'Contig'}->getSequence(0, $self->{'Exons'}->[-1]->[2], $self->{'Exons'}->[-1]->[2], 3, 0, $self->{'Strand'});
    $stop2 = substr($stop2, 1, 3);
#    print "$self->{'Id'} \"$stop1\" \"$stop2\"->", length($self->{'Seq'}), "\n";
    $self->{'NoStop'} = 1;
#    if(length($self->getFeatureSeq()) % 3 == 0 && 
    if($stop1 =~ /taa$/i || $stop1 =~ /tga$/i || $stop1 =~ /tag$/i) {
	$self->{'NoStop'} = 0;
	$self->trimStopCodon();
    }
#    if(length($self->getFeatureSeq()) % 3 == 0 && 
    elsif($stop2 =~ /taa$/i || $stop2 =~ /tga$/i || $stop2 =~ /tag$/i) {
	$self->{'NoStop'} = 0;
    } 
#    print "$self->{'NoStop'}\n";
    
}


sub trimStopCodon {
    my($self) = @_;
    my $lexon = $self->{'Exons'}->[-1];
    my$lenLE = abs($lexon->[2] - $lexon->[1]) + 1;

    if($self->{'Strand'} == STRAND_FWD) {
	if($lenLE > 3) {
	    $lexon->[2] -= 3;
	}
	else {
	    $lexon->[2] = $lexon->[1];
	    $self->{'NoStop'} = 1;
	}
	if(scalar(@{$self->{'3UTRExons'}}) > 0) {
	    $self->{'3UTRExons'}->[0]->[1] = $lexon->[2] + 1;
	}
    }
    else {
	if($lenLE > 3) {
	    $lexon->[2] += 3;
	}
	else {
	    $lexon->[2] = $lexon->[1];
	    $self->{'NoStop'} = 1;
	}
	if(scalar(@{$self->{'3UTRExons'}}) > 0) {
	    $self->{'3UTRExons'}->[0]->[1] = $lexon->[2] - 1;
	}
    }
    
    $self->getFeatureSeq();
}

sub computePhase {
    my($self) = @_;
    if(!defined($self->{'Phase'})) {
        $self->{'Phase'} = 0;
    }

    #find the first frame that contains no stop codons (or the least amount of them)
    my $phase;
    my $origSeq = $self->SUPER::getFeatureSeq();

    my $numStops = GenUtils::numPhaseOccurrences($origSeq, $self->{'Phase'}, "taa", "tga", "tag");
    if($numStops > 0) {
        print STDERR "$self->{'Contig'}->{'Id'} with $numStops stops in phase $self->{'Phase'}\t";
        ($phase, $numStops) = GenUtils::minPhaseOccurrences($origSeq, "taa", "tga", "tag");
        if(!$numStops) {
            $self->{'Phase'} = $phase;
            print STDERR "\tfixed with phase $phase\n";
            return 1;
        }
        print STDERR "\tnew phase $phase gives $numStops stops\n";
        return 0;

    }
    return 1;
}


sub getFeatureSeq {
    my($self) = @_;
    $self->{'Seq'} = $self->SUPER::getFeatureSeq();
    return $self->{'Seq'};
}

sub getThreePWeight {
    my $self = shift @_;
    return $self->getExonsWeight($self->{'3UTRExons'});
}

sub getFivePWeight {
    my $self = shift @_;
    return $self->getExonsWeight($self->{'5UTRExons'});
}

sub getCDSWeight {
    my $self = shift @_;
    return $self->getExonsWeight($self->{'Exons'});
}


sub displayHeader {
    my($self, $file, $contigs, $main_sep, $first_sep, $sec_sep) = @_;

    $self->displayId($file);
    my $fivePweights = $self->getWeights4Display($self->{'5UTRExons'});
    my $cdsweights = $self->getWeights4Display($self->{'Exons'});
    my $threePweights = $self->getWeights4Display($self->{'3UTRExons'});
    
    if(scalar(@{$self->{'5UTRExons'}}) > 0) {
        $self->displayLoc($file, $self->{'5UTRExons'}, $contigs,
			  $main_sep, $first_sep, $sec_sep);
        print $file "[";
    }

    if(defined($self->{'Id'}) && $self->{'NoStart'}) {
        if($self->{'Phase'} != 0) {
            print $file $self->{'Phase'};
        }
        print $file "<";
    }

    $self->displayLoc($file, $self->{'Exons'}, $contigs,
		      $main_sep, $first_sep, $sec_sep);
    if(defined($self->{'Id'}) && $self->{'NoStop'}) { 
        print $file ">";
    }

    if(scalar(@{$self->{'3UTRExons'}}) > 0) {
        print $file "]";
        $self->displayLoc($file, $self->{'3UTRExons'}, $contigs,
			  $main_sep, $first_sep, $sec_sep);
    }
    
    if($#{$threePweights} < 0 && $#{$fivePweights} < 0 &&
       $#{$cdsweights} < 0)  {
	return;
    }

    ($#{$fivePweights} == $#{$self->{'5UTRExons'}}) || die "Wrong number of 5 prime weights for $self->{'Id'}\n";
    ($#{$cdsweights} == $#{$self->{'Exons'}}) || die "Wrong number of CDS weights for $self->{'Id'}\n";
    ($#{$threePweights} == $#{$self->{'3UTRExons'}})  || die "Wrong number of 3 prime weights for $self->{'Id'}\n";

    print $file "\t";
    if($#{$fivePweights} >= 0)  {
	print $file join(";", @{$fivePweights});
	print $file "[";
    }
    if($#{$cdsweights} >= 0)  {
	print $file join(";", @{$cdsweights});
    }
    if($#{$threePweights} >= 0)  {
	print $file "]";
	print $file join(";", @{$threePweights});
    }
}

sub dump {
    my($self, $file, $contigs, $onlyHeader) = @_;
#    print "-->\"$self->{'Id'} $self->{'gId'} $self->{'tId'}\"\n";
    $onlyHeader = 0 if !defined($onlyHeader);
    $self->displayHeader($file, $contigs);

    print $file "\n";

    if(!$onlyHeader) {    
        $self->getFeatureSeq();
        $self->displaySeq($file);
    }
}

# TODO: make it work for genes with no coding regions, only UTRs
sub contains {
    my($self, $olapGene) = @_;
    my($exons_a, $exons_b) = ($self->{'Exons'}, $olapGene->{'Exons'});
    my $bool = $self->contigId() cmp $olapGene->contigId();
    my($min_a, $max_a, $min_b, $max_b) = (0,0,0,0);

    if($bool != 0 || $self->{'Strand'} != $olapGene->{'Strand'}) {
        return 0;
    }
    if(scalar(@$exons_a) == 0) {
        return 0;
    }
    if(scalar(@$exons_b) == 0) {
        return 1;
    }
    
    if($self->{'Strand'} == STRAND_FWD) {
        $min_a = $exons_a->[0]->[1];
        $max_a = $exons_a->[-1]->[2];
        $min_b = $exons_b->[0]->[1];
        $max_b = $exons_b->[-1]->[2];
    }
    else {
        $min_a = $exons_a->[-1]->[2];
        $max_a = $exons_a->[0]->[1];
        $min_b = $exons_b->[-1]->[2];
        $max_b = $exons_b->[0]->[1];
        
    }

    if($min_b < $min_a || $max_b > $max_a) {
        return 0;
    }
    
    my $TpFlank = -1;
    my $FpFlank = -1;
    my $i;
#    print STDERR "check ", $olapGene->{'NoStart'}, " ", $olapGene->{'NoStop'}, " ";
#    $olapGene->dump(\*STDERR, undef, 1);
#    print STDERR "against ", $self->{'NoStart'}, " ", $self->{'NoStop'}, " ";
#    $self->dump(\*STDERR, undef, 1);
    
    for($i = 0; $i < scalar(@$exons_a); $i++) {
        my $TpCond = ($exons_a->[$i]->[1] <= $exons_b->[0]->[1] && $exons_b->[0]->[1] <= $exons_a->[$i]->[2]);
        my $FpCond = ($exons_a->[$i]->[2] >= $exons_b->[-1]->[2] && $exons_b->[-1]->[2] >= $exons_a->[$i]->[1]);
        
        if($self->{'Strand'} == STRAND_COMP) {
            $TpCond = ($exons_a->[$i]->[1] >= $exons_b->[0]->[1] && $exons_b->[0]->[1] >= $exons_a->[$i]->[2]);
            $FpCond = ($exons_a->[$i]->[2] <= $exons_b->[-1]->[2] && $exons_b->[-1]->[2] <= $exons_a->[$i]->[1]);
        }
        
        if($TpFlank < 0 && $TpCond) {
            if(!$olapGene->{'NoStart'} && ($i != 0 || $exons_a->[$i]->[1] != $exons_b->[0]->[1])) {
                return 0;
            }
            $TpFlank = $i;
        }
        if($FpFlank < 0 && $FpCond) {
            if(!$olapGene->{'NoStop'} && (($i != scalar(@$exons_a) - 1) || $exons_a->[$i]->[2] != $exons_b->[-1]->[2])) {
                return 0;
            }
            $FpFlank = $i;
        }
    }
    if($TpFlank < 0 || $FpFlank < 0) {
        return 0;
    }
#    print STDERR "Flanks ", $TpFlank, " ", $FpFlank,"\n";
    for($i = $TpFlank; $i <= $FpFlank; $i++) {
#        print STDERR "\n $i $exons_a->[$i]->[1] != $exons_b->[$i - $TpFlank]->[1] $exons_a->[$i]->[2] != $exons_b->[$i - $TpFlank]->[2]";
        if($i != $TpFlank) {
            if($exons_a->[$i]->[1] != $exons_b->[$i - $TpFlank]->[1]) {
                return 0;
            }
        }
        if($i != $FpFlank) {
            if($exons_a->[$i]->[2] != $exons_b->[$i - $TpFlank]->[2]) {
                return 0;
            }
        }
    }
    return 1;
}

sub equalTo {
    my ($self, $feature) = @_;
    my $equal = 1;
    if(scalar(@{$self->{'5UTRExons'}}) && scalar(@{$feature->{'5UTRExons'}})) {
        $equal &&= $self->hasSameExonCoords($feature, '5UTRExons');
    }

    $equal &&=  $self->hasSameExonCoords($feature, 'Exons');
    
    if(scalar(@{$self->{'3UTRExons'}}) && scalar(@{$feature->{'3UTRExons'}})) {
        $equal &&= $self->hasSameExonCoords($feature, '3UTRExons');
    }
    return $equal;
}


sub shiftCoords {
    my($self, $shift, $contigLen, $newContigId) = @_;
    if(!defined($contigLen)) {
        $contigLen = $self->{'Contig'}->{'Len'};
    }
    if(!defined($newContigId)) {
        $newContigId = $self->contigId();
    }    

    if(scalar(@{$self->{'5UTRExons'}}) >= 0) {
        $self->shiftExonCoords($self->{'5UTRExons'}, $shift, $contigLen, $newContigId);
    }
    my @rmCode = @{$self->shiftExonCoords($self->{'Exons'}, $shift, $contigLen, $newContigId)};
    if($rmCode[1]) {
        $self->{'NoStart'} = 1;
        $self->{'Phase'} = ($self->{'Phase'} + $rmCode[1]) % 3;
    }
    if($rmCode[2]) {
        $self->{'NoStop'} = 1;
    }

    if(scalar(@{$self->{'3UTRExons'}}) >= 0) {
        $self->shiftExonCoords($self->{'3UTRExons'}, $shift, $contigLen, $newContigId);
    }
}

sub moveLocation2OrigContigs {
    my ($self, $contigs) = @_;
    my $exons = $self->{'Exons'};
    my $subcontig = $self->contigId();
#    print $subcontig, "\n";
    (defined($contigs->{$subcontig})) || die $self->{'Id'}." ".$subcontig;
    my $contig = $contigs->{$subcontig}->{'Exons'}->[0];
    $self->shiftCoords($contig->[1] - 1, 10000000000000, $contig->[0]);
}


sub cropCoords {
    my($self, $lowerL, $upperL) = @_;

    if(scalar(@{$self->{'5UTRExons'}}) >= 0) {
        $self->cropExonCoords($self->{'5UTRExons'}, $lowerL, $upperL);
    }

    my @rmCode = @{$self->cropExonCoords($self->{'Exons'}, $lowerL, $upperL)};
    if($rmCode[1]) {
        $self->{'NoStart'} = 1;
        $self->{'Phase'} = ($self->{'Phase'} + $rmCode[1]) % 3;
    }
    if($rmCode[2]) {
        $self->{'NoStop'} = 1;
    }
    if(scalar(@{$self->{'3UTRExons'}}) >= 0) {
        $self->cropExonCoords($self->{'3UTRExons'}, $lowerL, $upperL);
    }
#    print STDERR "left with ", scalar(@{$self->{'Exons'}}), " coding exons\n";
}


sub lowestCoord {
    my $self = shift @_;
    my $min = -1;
    if($self->{'Strand'} == STRAND_FWD) {
        if(scalar(@{$self->{'Exons'}})) {
            $min = $self->{'Exons'}->[0]->[1];
        }
        if(scalar(@{$self->{'5UTRExons'}}) > 0) {
            $min = $self->{'5UTRExons'}->[0]->[1];
        }
        elsif($min < 0 && scalar(@{$self->{'3UTRExons'}}) > 0) {
            $min = $self->{'3UTRExons'}->[0]->[1];
        }
	else {  return $min; }
    }
    else {
        if(scalar(@{$self->{'Exons'}})) {
            $min = $self->{'Exons'}->[-1]->[2];
        }
        if(scalar(@{$self->{'3UTRExons'}}) > 0) {
            $min = $self->{'3UTRExons'}->[-1]->[2];
        }
        elsif($min < 0 && scalar(@{$self->{'5UTRExons'}}) > 0) {
            $min = $self->{'5UTRExons'}->[-1]->[2];
	}
	else {  return $min; }
    }
    return $min;
}

sub highestCoord {
    my $self = shift @_;
    my $max = -1;
    if($self->{'Strand'} == STRAND_FWD) {
        if(scalar(@{$self->{'Exons'}})) {
            $max = $self->{'Exons'}->[-1]->[2];
#	    print STDERR "$self->{'Id'} max $max \"", join("_", @{$self->{'Exons'}->[0]}), "\" ", scalar(@{$self->{'Exons'}}), "\n";
        }
        if(scalar(@{$self->{'3UTRExons'}}) > 0) {
            $max = $self->{'3UTRExons'}->[-1]->[2];
        }      
        elsif($max < 0 && scalar(@{$self->{'5UTRExons'}}) > 0) {
            $max = $self->{'5UTRExons'}->[-1]->[2];
        }
	else { return $max; }
    }
    else {
        if(scalar(@{$self->{'Exons'}})) {
            $max = $self->{'Exons'}->[0]->[1];
        }
        if(scalar(@{$self->{'5UTRExons'}}) > 0) {
            $max = $self->{'5UTRExons'}->[0]->[1];
        }
        elsif($max < 0 && scalar(@{$self->{'3UTRExons'}})) {
            $max = $self->{'3UTRExons'}->[0]->[1];
        }
	else {  return $max; }
    }
    return $max;
}

sub contigId {
    my $self = shift @_;

    if(!scalar(@{$self->{'Exons'}})) {
        return scalar(@{$self->{'5UTRExons'}}) ? 
            $self->{'5UTRExons'}->[0]->[0]:
            $self->{'3UTRExons'}->[0]->[0];
    }
    
    return  $self->{'Exons'}->[0]->[0];
}

sub compSizeTo {
    my($a, $b) = @_;
    my ($contig_a, $contig_b);

    
    my $bool = $a->contigId() cmp $b->contigId();
    
    if($bool != 0) {
        return $bool;
    }

    my $len_a = $a->highestCoord() - $a->lowestCoord() + 1;
    my $len_b = $b->highestCoord() - $b->lowestCoord() + 1;

    return $len_a < $len_b;
}

sub fixAnnot4Alignment {
    my($self) = @_;
    my $exons = $self->{'Exons'};
    if(scalar(@{$exons}) == 0) {
	return;
    }

    my $strand = $self->{'Strand'};
    my $phase = $self->{'Phase'};
    
    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
	my $e_i = $exons->[$i];
	my $this_seq = $self->{'Contig'}->getSequence(0, $e_i->[1], $e_i->[2], 0, 1, $strand);
    
	my $numStops = GenUtils::numPhaseOccurrences($this_seq, $phase, "taa", "tga", "tag");
	if($numStops > 0) {
	    print STDERR "$self->{'Contig'}->{'Id'} $self->{'Id'} $e_i->[1] $e_i->[2] with $numStops stops in phase $phase\t";
	    ($phase, $numStops) = GenUtils::minPhaseOccurrences($this_seq, "taa", "tga", "tag");
	    if(!$numStops) {  print STDERR "\tfixed with phase $phase\n"; }
	    else { 
		print STDERR "\tnew phase $phase gives $numStops stops\n";
	    }
	}

	push(@{$exons->[$i]}, $phase);
	push(@{$exons->[$i]}, $numStops);

	my $seqExon = $self->{'Contig'}->getSequence(2, $e_i->[1], $e_i->[1], 2, 1, $strand);	
	if($i == 0 && !$self->{'NoStart'} && 
	   length($seqExon) >= 3 && $seqExon =~ /atg$/i) {
	    push(@{$exons->[$i]}, 1);
	}
	elsif($seqExon =~ /^ag/i) {
	    push(@{$exons->[$i]}, 1);
	}	
	else {
	    push(@{$exons->[$i]}, 0);	    
	}
	
	$seqExon = $self->{'Contig'}->getSequence(0, $e_i->[2], $e_i->[2], 3, 1, $strand);
	if($i == scalar(@$exons) - 1 && !$self->{'NoStop'} && length($seqExon) == 4 && 
	   ($seqExon =~ /\Staa/i || $seqExon =~ /\Stga/i || $seqExon =~ /\Stag/i)) {
	    ((($phase + abs($e_i->[2] - $e_i->[1]) + 1) % 3) == 0)
		|| warn "gene $self->{'Id'} $self->{'NoStop'} has bad stop phase\n";
	    push(@{$exons->[$i]}, 1);
	}
	elsif($seqExon =~ /^.gt/i) {
#($donor =~ /gt/i || $donor =~ /gc/i) || warn "$self->{'Id'} donor $seqExon $e_i->[2]\n";
	    push(@{$exons->[$i]}, 1);
	}
	else {
	    push(@{$exons->[$i]}, 0);	    
	}
	
	$phase = ($phase + abs($e_i->[2] - $e_i->[1]) + 1) % 3;
    }
}

sub compAlignQualityTo {
    my($a, $b) = @_;
    my ($contig_a, $contig_b);

    my $i;
    my $bool = $a->contigId() cmp $b->contigId();
    
    if($bool != 0) {
        return $bool;
    }
    my $trans_a = !$a->{'NoStart'} + !$a->{'NoStop'};
    my $trans_b = !$b->{'NoStart'} + !$b->{'NoStop'};

    if($trans_a != $trans_b) {
	return -1*($trans_a <=> $trans_b);
    }

    my $exons_a = $a->{'Exons'};
    my $exons_b = $b->{'Exons'};

    my $len_a = 0;
    my $numcan_a = 0;
    my $numstops_a = 0;
    my $avgq_a = 0;

    for($i = 0; $i < scalar(@{$exons_a}); $i++) {
	my $e = $exons_a->[$i];
	$len_a += abs($e->[2] - $e->[1]) + 1;
	$avgq_a += $e->[3];
	$numstops_a += $e->[5];
	$numcan_a += $e->[6] + $e->[7];
    }
    $numcan_a = 100*$numcan_a/($#$exons_a + 1);
    $numstops_a = 100*$numstops_a/($#$exons_a + 1);

    my $len_b = 0;
    my $numcan_b = 0;
    my $numstops_b = 0;
    my $avgq_b = 0;

    for($i = 0; $i < scalar(@{$exons_b}); $i++) {
	my $e = $exons_b->[$i];
	$len_b += abs($e->[2] - $e->[1]) + 1;
	$avgq_b += $e->[3];
	$numstops_b += $e->[5];
	$numcan_b += $e->[6] + $e->[7];
    }

    $numcan_b = 100*$numcan_b/($#$exons_b + 1);
    $numstops_b = 100*$numstops_b/($#$exons_b + 1);

#    if($numstops_a != $numstops_b) {
#	return $numstops_a <=> $numstops_b;
 #   }

#    if($numcan_a != $numcan_b) {
#	return -1*($numcan_a <=> $numcan_b);
#    }
        
#    if(int($avgq_a/20.0) != int($avgq_b/20.0)) {
#	return-1*( $avgq_a <=> $avgq_b);
#    }

    return  -1*($len_a + $numcan_a + $numstops_a + $avgq_a <=> $len_b + $numcan_b + $numstops_b + $avgq_b);
    
}

sub nexons_compareTo {
    my($a, $b) = @_;
    my $num_exons_a = scalar(@{$a->{'Exons'}});
    my $num_exons_b = scalar(@{$b->{'Exons'}});

    return $num_exons_a <=> $num_exons_b;

}

sub compareTo {
    my($a, $b) = @_;
    my ($contig_a, $contig_b);

    
    my $bool = $a->contigId() cmp $b->contigId();
    
    if($bool != 0) {
        return $bool;
    }

    my $i = 0;
    my $min_a = $a->lowestCoord();
    my $min_b = $b->lowestCoord();

    if($min_a != $min_b) {
        return $min_a <=> $min_b;
    }
    
    my $max_a = $a->highestCoord();
    my $max_b = $b->highestCoord();
    
    if($max_a != $max_b) {
	return $max_a <=> $max_b;
    }

    my $num_exons_a = $#{$a->{'Exons'}};
    my $num_exons_b = $#{$b->{'Exons'}};
    
    if($num_exons_a != $num_exons_b) {
	return $num_exons_a <=> $num_exons_b;
    }

    return $a->{'Strand'} <=> $b->{'Strand'};

}

sub addUTRAnnot {
    my($self, $mrna) = @_;
    my $has_utrs = 0;
    ($self->{'Strand'} == $mrna->{'Strand'}) ||
	die "strands for utr and loc $mrna->{'Id'} $self->{'Id'} are different\n";
    my $loc_exons = $self->{'Exons'};
    my $utr_exons = $mrna->{'Exons'};    
    
    my $fiveUTR = $self->{'5UTRExons'};
    my $j;
    for($j = 0; $j <= $#$utr_exons; $j++) {
	my $utr_e = $utr_exons->[$j];
	my $loc_e = $loc_exons->[0];
	if($self->{'Strand'} == STRAND_FWD) {
	    if($utr_e->[1] < $loc_e->[1]) {
		if($utr_e->[2] > $loc_e->[1]) {
		    push(@$fiveUTR, [$utr_e->[0], $utr_e->[1], $loc_e->[1] - 1]);
				     
		}
		else {
		    push(@$fiveUTR, $utr_e);		    
		}
		$has_utrs = 1;
	    }
	}
	if($self->{'Strand'} == STRAND_COMP) {
	    if($utr_e->[1] > $loc_e->[1]) {
		if($utr_e->[2] < $loc_e->[1]) {
		    push(@$fiveUTR, [$utr_e->[0], $utr_e->[1], $loc_e->[1] + 1]);
		}
		else {
		    push(@$fiveUTR, $utr_exons->[$j]);
		}
		$has_utrs = 1;
	    }
	}
    }
    
    my $threeUTR = $self->{'3UTRExons'};
    for($j = $#$utr_exons; $j >= 0; $j--) {
	my $utr_e = $utr_exons->[$j];
	my $loc_e = $loc_exons->[-1];
	if($self->{'Strand'} == STRAND_FWD) {
	    if($utr_e->[2] > $loc_e->[2]) {
		if($utr_e->[1] < $loc_e->[2]) {
		    unshift(@$threeUTR, [$utr_e->[0], $loc_e->[2]+1, $utr_e->[2]]);
		}
		else {
		    unshift(@$threeUTR, $utr_e);
		}
		$has_utrs = 1;
	    }
	}
	if($self->{'Strand'} == STRAND_COMP) {
	    if($utr_e->[2] < $loc_e->[2]) {
		if($utr_e->[1] > $loc_e->[2]) {
		    unshift(@$threeUTR, [$utr_e->[0], $loc_e->[2]-1, $utr_e->[2]]);
		}
		else {
		    unshift(@$threeUTR, $utr_e);
		}
		$has_utrs = 1;
	    }
	}

    }    

    return $has_utrs;
}

sub dumpUTRSeqs {
    my($self, $file) = @_;
    my $seq;
    my $fiveutr = $self->{'5UTRExons'};
    if(scalar(@{$fiveutr}) > 0) {
        print ">5UTR_", $self->{'Id'}, "\n";
        my $seq = $self->{'Contig'}->getSequence(0, $fiveutr->[0]->[1], $fiveutr->[-1]->[2], 0, 1, $self->{'Strand'});
        print $seq,"\n";
    }
    my $threeutr = $self->{'3UTRExons'};    
    if(scalar(@{$threeutr}) > 0) {
        print ">3UTR_", $self->{'Id'}, "\n";
        my $seq = $self->{'Contig'}->getSequence(0, $threeutr->[0]->[1], $threeutr->[-1]->[2], 0, 1, $self->{'Strand'});
        print $seq, "\n";
    }
}

1;
