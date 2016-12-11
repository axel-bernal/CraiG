package GenUtils;
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Utils;
use Gene;
use Organism;

require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
             computeGeneOverlaps
             computeGeneContainment
             loadSubContigLocs
             loadSubContigSeqs
             loadExtendedFasta
             createFilterValArrays
             printScoreFilterVals
             printSignalFilterVals
             loadFilterVals
             loadFilterVals4Id
             fasta2xfasta
             revComp
             readPhatGFFFormat
             readGFFFormat
             readLocsFormat
             readLocsFromString
             readContigLocation
             readUCSCFormat
             readGenBankFormat
             readGenScanFormat
             readGlimmerFormat
             writeGlimmerFormat
             writeGTFFormat
             displayContigs     
             displayGeneList   
             trimStopCodon
             extractIds
             revCompLocs
             parseGCClasses
             scramble
             split
             markIGenics
             extractFlankingIGenics
             produceAltSplices
             numPhaseOccurrences
             );

use strict;
### enumerates ###
use enum qw(STRAND_FWD=0 STRAND_COMP=1);


sub computeGeneOverlaps {
    my($genes, $rgenes, $olapGenes, $strandDep, $minOlap, $minWeight) = @_;
    my $gene;
    my $rgene;
    my $i = 0;
    my $j = 0;
#    print "hola ", scalar(@$genes), " ", scalar(@$rgenes), "\n";
    while($i < scalar(@$genes) && $j < scalar(@$rgenes)) {
	$rgene = $rgenes->[$j];
	$gene = $genes->[$i];
	my $r_beg = $rgene->lowestCoord();
	my $r_end = $rgene->highestCoord();
	my $beg = $gene->lowestCoord();
	my $end = $gene->highestCoord();
	my $r_strand = $rgene->{'Strand'};
	my $strand = $gene->{'Strand'};
	my $r_cid = $rgene->contigId();
	my $cid = $gene->contigId();
#	print STDERR "testing $gene->{'Id'} (", scalar(@{$gene->{'Exons'}}), " ", scalar(@{$rgene->{'Exons'}}), " ", $gene->lowestCoord(), "-", $gene->highestCoord(), " ", $cid, ") and $rgene->{'Id'} (", $rgene->lowestCoord(), "-", $rgene->highestCoord(), " ", $r_cid, ") \n";
	if($cid eq $r_cid) {
	    if($r_end - $beg < $minOlap || $minOlap > $end - $r_beg) {
		$i++ if $minOlap > $end - $r_beg;
		$j++ if $r_end - $beg < $minOlap;
		next;
	    }
	    else {
		my $r_weight = $rgene->{'Exons'}->[0]->[3] if scalar(@{$rgene->{'Exons'}});
		if((!defined $r_weight || $r_weight >= $minWeight) &&
		   ($strand == $r_strand || $strandDep != 1) &&
		   ($strand != $r_strand || $strandDep != -1)) {
#		    print "overlap $gene->{'Id'} $rgene->{'Id'} $rgene->{'Exons'}->[0]->[2] $i\n"; 
#		    push(@$olapGenes, splice(@$genes, $i, 1));
		    $olapGenes->{$i}++;
		    $i++;
		    next;
		}
		else {
		    $i++ if $r_end > $end;
		    $j++ if $end >= $r_end;
		    next;
		}
	    }
	}
	else {
	    my $b = $rgene->compareTo($gene);
	    $j++ if $b < 0;
	    $i++ if $b > 0;
	}
    }
}
	

sub computeGeneContainment {
    my($genes, $rgenes, $containment, $min_weight, $report_rlocs) = @_;
    my $gene;
    my $rgene;
    my $i = 0;
    my $j = 0;
#    print "hola ", scalar(@$genes), " ", scalar(@$rgenes), "\n";
    while($i < scalar(@$genes) && $j < scalar(@$rgenes)) {
	$rgene = $rgenes->[$j];
	$gene = $genes->[$i];
	my $r_beg = $rgene->lowestCoord();
	my $r_end = $rgene->highestCoord();
	my $beg = $gene->lowestCoord();
	my $end = $gene->highestCoord();
	my $r_strand = $rgene->{'Strand'};
	my $strand = $gene->{'Strand'};
	my $r_cid = $rgene->contigId();
	my $cid = $gene->contigId();
#	print "testing $gene->{'Id'} ($cid"."_", $beg, "-", $end,  ") and $rgene->{'Id'} ($r_cid"."_", $r_beg, "-", $r_end, ")\n";
	if($strand == $r_strand) {
	    if($r_end < $beg || $r_beg > $end) {
		$i++ if $r_beg > $end;
		$j++ if $r_end < $beg;
		next;
	    }
	    else {
		my $r_exons = $rgene->{'Exons'};
		my $exons = $gene->{'Exons'};
		my $r_weight = $r_exons->[0]->[3] if scalar(@{$r_exons});

		if(defined $r_weight && $r_weight < $min_weight) {
		    $j++;
		    next;
		}

		if($#$exons == 0 || $#$r_exons == 0) {
		    $i++ if $#$exons == 0;
		    $j++ if $#$r_exons == 0;
		    next;
		}

#		print "overlap\n";
		my %hr_exons = ();
		for(my $i = 0; $i <= $#$r_exons; $i++) {
		    if($i) {
			my $loc_i = $r_exons->[$i - 1]->[2]."_".$r_exons->[$i]->[1];
			$hr_exons{$loc_i} = $i;
		    }
		}

		my $last_exon = -1;
		my $contained = 1;
		for(my $i = 0; $i <= $#$exons; $i++) {
		    if($i) {
			my $loc_i = $exons->[$i - 1]->[2]."_".$exons->[$i]->[1];
			if(!defined($hr_exons{$loc_i})) {
			    $contained = 0;
			    last;
			}
#			if($last_exon != -1 && $hr_exons{$loc_i} != $last_exon + 1) {
#			    $contained = 0;
#			    last;
#			}
			$last_exon = $hr_exons{$loc_i};
		    }
		}

		if($contained) {
#		    print "$gene->{'Id'} is contained in $rgene->{'Id'} $rgene->{'Exons'}->[0]->[2] $i\n"; 
#		    push(@$olapGenes, splice(@$genes, $i, 1));
		    if($report_rlocs) {
			$containment->{$j}++;
			$j++;			
		    }
		    else {
			$containment->{$i}++;
			$i++;
		    }
		    next;
		}
		else {
		    if($report_rlocs) {
			$j++;
			$i = 0;
		    }
		    else {
			$j = 0;
			$i++;
		    }
		    next;
		}
	    }
	}
	else {
	    my $b = $rgene->compareTo($gene);
	    $j++ if $b < 0;
	    $i++ if $b > 0;
	}
    }
}
	

sub loadSubContigLocs {
    my $file = shift @_;
    my $line;
    my %contigRanges;
    my %subcontigLengths;

    while(defined($line = <$file>) && ($line =~ /^>(\S+)\s+(.+)\s+?(\d+)?b?p?$/)) {
        my $id = $1;
        my @locs = @{&GenUtils::readLocExons($2)};
        my $len = 0;
        foreach my $exon (@locs) {
            my $cid = $exon->[0];
            if(!defined($contigRanges{$cid})) {
                @{$contigRanges{$cid}} = ();
            }

            my $i = 0;            
#            print "\t\t$cid $id [ $len $exon->[1] $exon->[2] ", 
#            scalar(@{$contigRanges{$cid}}), "\t";
            
            while($i < scalar(@{$contigRanges{$cid}}) &&
                  $contigRanges{$cid}->[$i]->[2] < $exon->[1]) {
                $i++;
            }
            splice(@{$contigRanges{$cid}}, $i, 0, [$id, $len, $exon->[1], $exon->[2]]);

#            print scalar(@{$contigRanges{$cid}}), "\n";
            $len += abs($exon->[2] - $exon->[1]) + 1;
        }

        $subcontigLengths{$id} = $len;
    }    

    return [\%contigRanges, \%subcontigLengths];
}


sub loadSubContigSeqs {
    my $file = shift @_;
    my $org = Organism->new();
    my %subContigs = % {$org->loadContigs(\*CONTIGS)};
    my %subcontigLengths = ();
    my %contigRanges = ();

    foreach my $id (keys %subContigs) {
        my $offset = 1;
        my $cid = $id;
        if($id =~ /^(\S+)\.(\d+)$/) {
            $offset = $2;
            $cid = $1;
        }

        my $len = $subContigs{$id}->{'Len'};
        
        if(!defined($contigRanges{$cid})) {
            @{$contigRanges{$cid}} = ();
        }
        
        $subcontigLengths{$id} = $len;

        my $i = 0;            
        while($i < scalar(@{$contigRanges{$cid}}) &&
              $contigRanges{$cid}->[$i]->[2] < $offset) {
            $i++;
        }
        splice(@{$contigRanges{$cid}}, $i, 0, [$id, 0, $offset, $offset + $len - 1]);
    }
    return [\%contigRanges, \%subcontigLengths];

}
sub loadExtendedFasta {
    my $file = shift @_;
    my $line = <$file>;
    my %marks = ();
    
    while(defined($line) && $line =~ /^>(\S+)\s*/) {
        my $id = $1;
        my @seq = ();
        while(defined($line = <$file>) && $line !~ /^>/) {
            $line =~ /^([^\s,>])\s+x\s+(\d+)/;
            push(@seq, [$1, $2]);
        }
        $marks{$id} = \@seq;
    }
    return \%marks;
}


sub fasta2xfasta {
    my ($id, $seq_ref) = @_;
    my $result = "";
    if(length ($seq_ref) < 1) {
        return;
    }

    my $offset = 1;
    if($id =~ /^\S+\.(\d+)$/) {
        $offset = $1;
    }

    my $last_char = substr($seq_ref, 0, 1);
    my $times = 1;
    my $pos;
    for(my $i=1;$i< length($seq_ref);$i++)
    {	
        my $this_char = substr($seq_ref, $i, 1);
        if($this_char ne $last_char) {
            $pos =  $offset + $i - $times;
            $result .= $last_char." x $times $pos\n";
            $last_char = $this_char;
            $times = 1;
        }
        else {
            $times++;
        }
        
    }
    $pos = $offset + length($seq_ref) - $times;
    $result .= $last_char." x $times $pos\n";
    return $result;
}

sub createFilterValArrays  {
    my($contigScores, $c_id, $length, $def_val) = @_;
    my $scores;
    if(!defined($def_val)) {
	$def_val = 0;
    }
    if(!defined($contigScores->{$c_id})) {
	my @scores_fwd = ();
	my @scores_bwd = ();
	$contigScores->{$c_id} = [\@scores_fwd, \@scores_bwd];
	$scores = $contigScores->{$c_id};
	
	for(my $j = 0; $j <= $length; $j++) {
	    push(@{$scores->[0]}, $def_val);
	    push(@{$scores->[1]}, $def_val);
	}
    }
    else {
	$scores = $contigScores->{$c_id};
    }

    return $scores;
}


sub loadFilterVals {
    my($file, $def_val) = @_;
    my $line =  <$file>;
    my %filterVals;
    my $cid = "";
    while(1) {
	$cid = loadFilterVals4Id($file, \$line, \%filterVals, $def_val);
	if($cid eq "") {
	    last;
	}
    }
    return \%filterVals;
}

sub loadFilterVals4Id {
    my($file,$ref_line, $filterVals, $def_val) = @_;
    my $cid = "";
    if(defined($$ref_line) && $$ref_line  =~ />(\S+)\s+(\d+)/) {
	$cid = $1;
	my $clen = $2;
	my $strand = STRAND_FWD;
#	print "$cid\n";
	my $fVals = GenUtils::createFilterValArrays($filterVals, $cid, $clen, $def_val);

	my $sigPos1 = 1;
	
	while(defined($$ref_line = <$file>) && $$ref_line  !~ />(\S+)\s+(\d+)/) {
	    my $sigPos2 = 1;
	    my $sigVal = "";
	    if($$ref_line =~ /^(\S)$/) { # xfasta
		$sigVal = $1;
		$filterVals->{$cid}->[$strand]->[$sigPos1++] = $sigVal;
	    }
	    if($$ref_line =~ /^(\S)\s+x\s+(\d+)(\s*.*)$/) { #signals or xfasta
		$sigPos2 = $sigPos1 + $2 - 1;
		$sigVal = "$1$3";
		chomp $sigVal;
		for(my $sigPos = $sigPos1; $sigPos <= $sigPos2; $sigPos++) {
		    $filterVals->{$cid}->[$strand]->[$sigPos] = $sigVal;
		}
#		print "$strand $sigPos1..$sigPos2 $sigVal $fVals->[$strand]->[$sigPos1]\n";
		$sigPos1 = $sigPos2 + 1;
	    }
	    elsif($$ref_line =~ /^([\d,\.]+)\s+([\-,\d,.]+)\s*$/) { #scores
		$sigPos1 = $sigPos2 = $1;
		$sigVal = $2;
		if($1 =~ /^([^\.]+)\.\.(\d+)$/) {
		    $sigPos1 = $1;
		    $sigPos2 = $2;
		}
		for(my $sigPos = $sigPos1; $sigPos <= $sigPos2; $sigPos++) {
#		    print STDERR "$sigPos $sigVal\n" if $sigPos == 415;
		    $fVals->[$strand]->[$sigPos] += $sigVal;                
		}
	    }
	    elsif($$ref_line =~ /\/\//) {
		$sigPos1 = 1;
		$strand = STRAND_COMP;
	    }
	}
	
	if($strand == STRAND_FWD) {
	    # this means that the filter did not contain a comp strand info and so the
	    # comp strand has reversed info from the fwd.
	    for(my $sigPos = 1; $sigPos <= $clen; $sigPos++) {
		$fVals->[STRAND_COMP]->[$sigPos] = $fVals->[$strand]->[$clen - $sigPos + 1];
	    }
	}
    }
    return $cid;
}


sub printContigScoreFilterVals {
    my($cid, $filter_vals, $file, $strand) = @_;
    my $fvals = $filter_vals->{$cid};
    my $empty = 1;

    for(my $t = 0; $t < 2; $t++) {
	my $length = scalar(@{$fvals->[$t]});

	if(!$t) {
	    print $file ">$cid\t", $length - 1, "\n";
	}
	
	if($t && ($strand == STRAND_COMP || $empty)) {
	    print $file "//\n";
	}
	
	my $i = 2;
	my $last_fval = $fvals->[$t]->[1];
	my $last_pos = 1;
	
	for( ; $i < $length; $i++) {
	    if($fvals->[$t]->[$i] ne $last_fval) {
		if($last_fval != 0) {
		    print $file $last_pos;
		    if($last_pos < $i - 1) {
			print $file "..", $i - 1;
		    }
		    print $file "\t", $last_fval, "\n";
		    $empty = 0;
		}
		$last_fval = $fvals->[$t]->[$i];
		$last_pos = $i;
	    }
	}
#       print STDERR "\"$t $cid \"$last_fval\" $last_pos \"\n";
	if($last_fval != 0) {
	    print $file $last_pos;
	    if($last_pos < $length - 1) {
		print $file "..", $length - 1;
	    }
	    print $file "\t", $last_fval, "\n";
	    $empty = 0;
	}
    }
}


sub printScoreFilterVals {
    my ($filter_vals, $file, $strand) = @_;

    foreach my $cid (keys %$filter_vals) {
	printContigScoreFilterVals($cid, $filter_vals, $file, $strand);
    }
}

sub printSignalFilterVals {
    my ($filter_vals, $file, $strand) = @_;
    
    foreach my $cid (keys %$filter_vals) {
        for(my $t = 0; $t < 2; $t++) {
	    my $fvals = $filter_vals->{$cid}->[$t];
            my $length = scalar(@{$fvals});

            if(!$t) {
                print $file ">$cid\t", $length - 1, "\t3\n";
            }

            if($t && $strand == STRAND_COMP) {
                print $file "//\n";
            }

            my $i = 2;
	    my $last_fval = $fvals->[1];
	    my $last_pos = 1;
            for( ; $i < $length; $i++) {
                if($fvals->[$i] ne $last_fval) {
		    $last_fval =~ /^(\S)(\s*.*)$/;
		    print $file $1, " x ", $i - $last_pos, 
		    (defined($2) && length($2) > 0 ? $2 :""), "\n";
#		    print "hola $1 $2 \"$last_fval\"\n";

		    $last_fval = $fvals->[$i];
		    $last_pos = $i;
                }
            }

	    if($last_pos < $length) {
		$last_fval =~ /^(\S)(\s*.*)$/;
		print $file $1, " x ", $length - $last_pos, 
		(defined($2) && length($2) > 0 ? $2 :""), "\n";
#		print "hola $1 $2 \"$last_fval\"\n";
	    }
        }
    }
}

sub fillSigFastaArray {
    my($gene, $weight, $array, $symbol, $pos, $key, $strand) = @_;
#    print STDERR "new sig $symbol @ $pos $key $weight\n";
    my $parray = $array->[$strand][$pos];
    my $conflict = 0;
    if($parray->[0] ne "-") {
        if($parray->[0] ne $symbol) {
	    $conflict = 1;
	    warn "$gene->{'Id'} \@ position $pos has two different signals\n";
	    my $oldw = 0;
	    for(my $k = 1; $k < scalar(@{$parray}); $k++) {
		$oldw = $parray->[$k]->[1] unless $oldw > $parray->[$k]->[1];
	    }
	    warn "existing weight=$oldw, new weight=$weight\n";
	    if(($symbol ne "S" && $symbol ne "T") || $oldw > $weight) {
		return 0;
	    }
	    warn "discarding the one with lower weight\n";
	}
    }
    
    if($conflict) {
	@{$parray} = ($symbol, [$key, $weight]) 
    }
    else {
	$parray->[0] = $symbol;
	my $found = 0;
	for(my $k = 1; $k < scalar(@{$parray}); $k++) {
	    if($key eq $parray->[$k]->[0]) {
		$parray->[$k]->[1] += $weight;
		$found = 1;
		last;
	    }
	}
	push(@{$parray}, [$key, $weight]) unless $found;
    }
    return 1;
}

sub fillSigScoreArray {
    my($gene, $weight, $countArray, $maxArray, $pos, $maxCounts, $strand) = @_;

    if($maxArray->[$strand][$pos] < $weight) {
	$maxArray->[$strand][$pos] = $weight;
    }

    $countArray->[$strand][$pos]++;    

    if($countArray->[$strand][$pos] > $$maxCounts) {
	$$maxCounts = $countArray->[$strand][$pos];
    }

}

sub fillSigArray {
    my($gene, $weight, $countArray, $maxArray, $array, $symbol, $pos, $offsetP, $cLen, $lastPos, $lastOffset, $phase, $maxCounts, $strand) = @_;
    my $cpos = (($strand == STRAND_FWD) ? $pos + $offsetP: $cLen - ($pos - $offsetP) + 1);
    if($lastPos > 0) {
        $lastPos = (($strand == STRAND_FWD) ? $lastPos + $lastOffset: $cLen - ($lastPos - $lastOffset) + 1);
    }

    my $key = "$lastPos $phase";

    if(fillSigFastaArray($gene, $weight, $array, $symbol, $cpos, $key, $strand)) {
        fillSigScoreArray($gene, $weight, $countArray, $maxArray, $cpos, $maxCounts, $strand);
	return 1;
    }
    return 0;
    
}


sub printSigArray {
    my ($a_ref, $file) = @_;
    my $k;

    if(scalar(@$a_ref) < 2) {
        return $a_ref->[1];
    }
    my $last_char = $a_ref->[1]->[0];
    my $last_content = "";
    for($k = 1; $k < scalar(@{$a_ref->[1]}); $k++) {
	$last_content .= join(" ", @{$a_ref->[1]->[$k]});
	$last_content .= " " unless $k == scalar(@{$a_ref->[1]}) - 1;
    }
    
    my $times = 1;
    for(my $i=2;$i< scalar @$a_ref; $i++) {
        my $this_char = $a_ref->[$i]->[0];
	my $this_content = "";
	for($k = 1; $k < scalar(@{$a_ref->[$i]}); $k++) {
	    $this_content .= join(" ", @{$a_ref->[$i]->[$k]});
	    $this_content .= " " unless $k == scalar(@{$a_ref->[$i]}) - 1;
	}

        if($this_char ne $last_char || $this_content ne $last_content) {
            print $file $last_char," x $times";
	    print $file $last_content ne "" ? " ".$last_content : "", "\n";
            $last_char = $this_char;
	    $last_content = $this_content;
            $times = 1;
        }
        else {
            $times++;
        }
        
    }

    print $file $last_char," x $times";
    print $file $last_content ne "" ? " ".$last_content : "", "\n";
}


sub printStateArray {
    my ($seq_ref, $out_stream) = @_;

    if(scalar @{$seq_ref} < 2) {
        return;
    }

    my $last_char = $seq_ref->[1];
    my $times = 1;
    my $pos; 
    my $need_lc = 0;

    for(my $i=2;$i<scalar @{$seq_ref};$i++)
    {	
        if($seq_ref->[$i] ne $last_char) {
            $pos = $i - $times;
            if($times > 10) {
                if($need_lc) {
                    print $out_stream "\n";
                    $need_lc = 0;
                }
                print $out_stream $last_char," x $times\t$pos\n";
            }
            else {
                for(my $j = 0; $j < $times; $j++)  {
                    print $out_stream $last_char;
                }
                $need_lc = 1;
            }
            $last_char = $seq_ref->[$i];
            $times = 1;
        }
        else {
            $times++;
        }
        
    }
    $pos = scalar(@{$seq_ref}) - $times;
    if($times > 10) {
        if($need_lc) {
            print $out_stream "\n";
        }
        print $out_stream $last_char, " x $times\t$pos\n";   
    }
    else {
        for(my $j = 0; $j < $times; $j++)  {
            print $out_stream $last_char;
        }
        print $out_stream "\n";
    }
}


sub revComp {
    my( $seqP ) = @_;
    my( $rev  );
    
    $rev =  reverse( $seqP );
    $rev =~ tr/ACGTUMRWSYKBDHVacgtumrwsykbdhv/TGCAAKYWSRMVHDBtgcaakywsrmvhdb/;
    return $rev;
}

#Format is AB000381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";

sub readPhatGFFFormat {
    my ($file, $contigs) = @_;
    my $line = "";
    my %genes;
    my ($geneId, $oldGeneId) = ("", "");
    #handling eukaryotic genes
    my @exons = ();
    my $phase = 0;
    my @entry;    
    while(defined($line = <$file>)) {
	@entry = ();
	my $exon;

        if($line =~ /^\#/) {
            next;
        }
        @entry = split(/\s+/,$line);
        # GFF goes like this : <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments], [attributes] usually have gene "gene_id" as first attribute
	
	$entry[9] =~ /\"?([^\"]+)\"?\;?/;
	$geneId = $entry[0]."_".$1;
	if($entry[2] ne "CDS"){
	    next;
	}
	my $begExon = $entry[3];
	my $endExon = $entry[4];
	if($entry[6] eq "-") {
	    $begExon = $contigs->{$geneId}->{'Len'} - $begExon + 1;
	    $endExon = $contigs->{$geneId}->{'Len'} - $endExon + 1;
	}
	$exon =  [$entry[0], $begExon, $endExon];
	if($oldGeneId ne "" && $geneId ne $oldGeneId) { # We have changed gene		    
	    my @tmp = @exons;
	    my $gene = Gene->new('Id' => $oldGeneId, 'Exons' => \@tmp, 'Seq' => "", 'Strand' => ($entry[6] eq "+" ? STRAND_FWD:STRAND_COMP), 'Phase' => $phase);
	    @exons = ();
            $phase = 0;
	    $genes{$oldGeneId} = $gene;
	}
        if(scalar(@exons) == 0) {
            $phase = $entry[7];
        }
	push(@exons, [$entry[0], $begExon, $endExon]);
	$oldGeneId = $geneId;
    }
    # Inserting last gene
    my $gene = Gene->new('Id' => $oldGeneId, 'Exons' => \@exons, 'Seq' => "", 'Strand' => ($entry[6] eq "+" ? STRAND_FWD:STRAND_COMP), 'Phase' => $phase);
    $genes{$geneId} = $gene;    
    return \%genes;
}

sub checkAnnotationSanity {
    my($genes, $preserve_phases) = @_;
    fixTranslationProduct($genes, $preserve_phases);
    checkTranscribedProduct($genes);
}

sub fixTranslatedProduct {
    my($genes, $preserve_phases) = @_;
#    print scalar keys %$genes, " \n";
    $preserve_phases = 0 if !defined($preserve_phases);
    foreach my $key (keys %$genes) {
	my $gene = $genes->{$key};
        $gene->checkStop();
	$gene->computePhase() if !$preserve_phases;
	$gene->checkStart();
    }

}

sub checkTranscribedProduct {
    my($genes) = @_;
    my $donor;
    my $acceptor;
    my $exon_i;
    my $exon_i_1;

#	print "INFO $gene->{'Id'} $gene->{'NoStart'} $gene->{'NoStop'}\n";
    foreach my $key (keys %$genes) {
	my $gene = $genes->{$key};
    
	for(my $i = 0; $i < scalar(@{$gene->{'Exons'}}); $i++) {
	    $exon_i = $gene->{'Exons'}->[$i];
	    if($i > 0) {
		my $exon_i_1 = $gene->{'Exons'}->[$i - 1];
		$donor = substr($gene->{'Contig'}->getSequence(0, $exon_i_1->[1], $exon_i_1->[2], 2, 1, $gene->{'Strand'}), abs($exon_i_1->[2] - $exon_i_1->[1]) + 1, 2);
		$acceptor = substr($gene->{'Contig'}->getSequence(2, $exon_i->[1], $exon_i->[2], 0, 1, $gene->{'Strand'}), 0, 2);                
		
		if(abs($exon_i->[1] - $exon_i_1->[2]) < 4) {
		    if($donor !~ /gt/i || $acceptor !~ /ag/i) {
			warn "Small intron (<4) $gene->{'Id'}, $i, $exon_i->[1], $exon_i_1->[2], with no good splice signals $donor $acceptor\n";
			$exon_i_1->[2] = $exon_i->[2];
			splice(@{$gene->{'Exons'}}, $i, 1);
			$i--;
			next;
		    }
		}
		
		if($donor !~ /gt/i || $acceptor !~ /ag/i) {
		    warn "Intron $gene->{'Id'}, $i, $exon_i->[1], $exon_i_1->[2], with no good splice signals $donor $acceptor\n";                    
		}
	    }
	}
    }
}


sub fixAnnotation {
    my($genes) = @_;
    my $donor;
    my $acceptor;
    my $exon_i;
    my $exon_i_1;

    foreach my $key (keys %$genes) {
        my $gene = $genes->{$key};
        my $state = 0;

        for(my $i = 0; $i < scalar(@{$gene->{'Exons'}}); $i++) {
            $exon_i = $gene->{'Exons'}->[$i];
            if($i > 0) {
                my $exon_i_1 = $gene->{'Exons'}->[$i - 1];
                $donor = $gene->{'Contig'}->getSequence(0, $exon_i_1->[2], $exon_i_1->[2], 3, 1, $gene->{'Strand'});
                $acceptor = $gene->{'Contig'}->getSequence(3, $exon_i->[1], $exon_i->[1], 0, 1, $gene->{'Strand'});
                
                if($donor !~ /.gt./i) {
                    if($donor =~ /gt../i) {
                        if($state < 0) {
                            $exon_i_1->[1] += $gene->{'Strand'} == STRAND_FWD ? -1 : 1;
                            $state = -2;
                        }
                        else {
                            $state = -1;
                        }
                    }
                    elsif($donor =~ /..gt/i) {
                        if($state > 0) {
                            $exon_i_1->[1] += $gene->{'Strand'} == STRAND_FWD ? 1 : -1;
                            $state = 2;
                        }
                        else {
                            $state = 1;
                        }
                    }
                }
                else {
                    if($state == 2) {
                        $exon_i_1->[1] += $gene->{'Strand'} == STRAND_FWD ? 1 : -1;
                        $state = 0;
                    }
                    elsif($state == -2) {
                        $exon_i_1->[1] += $gene->{'Strand'} == STRAND_FWD ? -1 : 1;
                        $state = 0;
                    }
                }

                if($acceptor !~ /.ag./i) {
                    if($acceptor =~ /ag../) {
                        if($state < 0) {
                            $exon_i_1->[2] += $gene->{'Strand'} == STRAND_FWD ? -1 : 1;
                            $state = -2;
                        }
                        else {
                            $state = -1;
                        }
                    }
                    elsif($acceptor =~ /..ag/) {
                        if($state > 0) {
                            $exon_i_1->[2] += $gene->{'Strand'} == STRAND_FWD ? 1 : -1;
                            $state = 2;
                        }
                        else {
                            $state = 1;
                        }                        
                    }                    
                }
                else {
                    if($state == 2) {
                        $exon_i_1->[2] += $gene->{'Strand'} == STRAND_FWD ? 1 : -1;
                        $state = 0;
                    }
                    elsif($state == -2) {
                        $exon_i_1->[2] += $gene->{'Strand'} == STRAND_FWD ? -1 : 1;
                        $state = 0;
                    }
                }
            }
        }
    }
}


sub extractUTRExons {
    my($utrExons, $exons, $strand) = @_;
     my @fiveUtrExons = ();
    my @threeUtrExons = ();
    my $firstExons = \@fiveUtrExons;
    my $lastExons = \@threeUtrExons;

    if($strand == STRAND_COMP) {
        $firstExons = \@threeUtrExons;
        $lastExons = \@fiveUtrExons;
    }

    for(my $i = 0; $i < scalar(@{$utrExons}); $i++) {
        if(scalar(@$exons) != 0) {
	   if($utrExons->[$i]->[2] > $exons->[0]->[1]) {
	       push(@$lastExons, $utrExons->[$i]);
	   }
	   else {
	       push(@$firstExons, $utrExons->[$i]);
	   }
	}
	else {
	    push(@$firstExons, $utrExons->[$i]);
	}
    }

    return [\@fiveUtrExons, \@threeUtrExons];
}

sub insertGFFGene {
    my($genes, $exons, $utrExons, $lengths, $strand, $noStart, $noStop, $phase, $onlyOneTranscript, $contig, $oldGeneId, $type, $addStart, $preserve_phases) = @_;
    
    $preserve_phases = 0 if !defined($preserve_phases);

    my $tId;
    if($onlyOneTranscript) {
        # choose the longest transcript or the one with longer coding regions if 
	# equal length
        my $maxCodingLength = 0;
        my $chosentId = undef;
        my $maxLength = 0;

        foreach $tId (keys %{$exons}) {
	    my $transcriptLength = abs($exons->{$tId}->[-1]->[2] - $exons->{$tId}->[0]->[1]) + 1;
	    if($transcriptLength > $maxLength) {
		$chosentId = $tId;
		$maxLength = $transcriptLength;
	    }
            elsif($transcriptLength == $maxLength) {
		if($lengths->{$tId} > $maxCodingLength) {
		    $chosentId = $tId;
		    $maxCodingLength = $lengths->{$tId};
		}
            }
        }

	if(defined($tId)) {
	    my $tmp_exons = $exons->{$tId};
	    %{$exons} = ();
	    $exons->{$tId} = $tmp_exons;
	}
    }

    foreach $tId (keys %{$exons}) {
        my $transId = $tId;
        if($oldGeneId ne $transId) {
            $transId = $oldGeneId."|".$tId;
        }

#        print "**\"$transId\" \"$oldGeneId\" $tId $contig ",scalar(@{$exons->{$tId}}),"\n";
        
        my $gene;

        my @tmp = @{$exons->{$tId}};
        if(scalar(@tmp)) {
#            print STDERR "inserting ", scalar(@tmp), " exons for gene $transId ", $noStart->{$tId}, " ", $noStop->{$tId}, "\n";

            if(!defined($genes->{$transId})) {
                $gene = Gene->new('tId' => $tId, 'gId' => $oldGeneId, 'Exons' => \@tmp, 'Contig' => $contig, 'Seq' => "", 'Strand' => $strand->{$tId}, 'NoStart' => $noStart->{$tId}, 'NoStop' => $noStop->{$tId}, 'Phase' => $phase->{$tId});
            }
            else { # transcript already inserted
                $gene = $genes->{$transId};
                
                for(my $i = 0; $i < scalar(@tmp); $i++) {
                    if($strand->{$tId} == STRAND_FWD) {
                        push(@{$genes->{$transId}->{'Exons'}}, $tmp[$i]);
                    }
                    else {
                        unshift(@{$genes->{$transId}->{'Exons'}}, $tmp[scalar(@tmp) - $i - 1]);
                        $genes->{$transId}->{'Phase'} = $phase->{$tId};
                    }
                }
                $gene->getFeatureSeq();
            }
            
            if(defined($contig)) {
                $contig->{'Features'}->{$type}->{$transId} = $gene;
		$gene->checkStop();
		$gene->computePhase() if !$preserve_phases;
		$gene->checkStart();
            }
            
            $genes->{$transId} = $gene;

        }
    }

    foreach $tId (keys %{$utrExons}) {
        my $transId = $tId;
        if($oldGeneId ne $transId) {
            $transId = $oldGeneId."|".$tId;
        }
#                print "**\"$transId\" \"$oldGeneId\" $tId $contigKey $entry[2] ",scalar(@{$exons{$tId}}),"\n";
        
        if(!scalar(@{$utrExons->{$tId}})) {
            next;
        }

        my $gene;
        
        if(!defined($genes->{$transId})) {
            my($tmp1, $tmp2) = @{&GenUtils::extractUTRExons($utrExons->{$tId}, $exons->{$tId}, $strand->{$tId})};
#            print STDERR "inserting ", scalar(@$tmp1), " ", scalar(@$tmp2), " UTR exons for new gene $transId\n";

            $gene = Gene->new('tId' => $tId, 'gId' => $oldGeneId, '5UTRExons' => \@$tmp1, '3UTRExons' => \@$tmp2, 'Contig' => $contig, 'Seq' => "", 'Strand' => $strand->{$tId}, 'NoStart' => $noStart->{$tId}, 'NoStop' => $noStop->{$tId}, 'Phase' => $phase->{$tId});
        }
        else { # transcript already inserted
            $gene = $genes->{$transId};
            my($tmp1, $tmp2) = @{&GenUtils::extractUTRExons($utrExons->{$tId}, $gene->{'Exons'}, $strand->{$tId})};
#            print "inserting ", scalar(@$tmp1), " ", scalar(@$tmp2), " UTR exons for existing gene $transId\n";

            for(my $i = 0; $i < scalar(@$tmp1); $i++) {
                if($strand->{$tId} == STRAND_FWD) {
                    push(@{$gene->{'5UTRExons'}}, $tmp1->[$i]);
                }
                else {
                    unshift(@{$gene->{'5UTRExons'}}, $tmp1->[scalar(@$tmp1) - $i - 1]);
                }
            }

            for(my $i = 0; $i < scalar(@$tmp2); $i++) {
                if($strand->{$tId} == STRAND_FWD) {
                    push(@{$gene->{'3UTRExons'}}, $tmp2->[$i]);
                }
                else {
                    unshift(@{$gene->{'3UTRExons'}}, $tmp2->[scalar(@$tmp2) - $i - 1]);
                }
            }
        }

        $genes->{$transId} = $gene;
    }
}

#
# Reads tigr and twinscan format. Commented code looks for contig id in 
# twinscan format.
#
sub readGFFFormat {
    my ($file, $type, $contigs, $onlyOneTranscript, $uniq_id, $preserve_phases) = @_;
    
    $preserve_phases = 0 if !defined($preserve_phases);
    my $line = "";
    my %genes;
    my $numExons = 0;
    my $geneId = "";
    my $oldGeneId = "";
    my $transcriptId = "";
    my %utrExons = ();
    my %exons = ();
    my @entry;    
    my %noStop = ();
    my %noStart = ();
    my %strand = ();
    my %lengths = ();
    my %phase = ();
    my $exons;
    my $contig;
    my $contigKey;
    my $counter = 0;
    my $addStart = 1;
    my $source = "";

    if(!defined($onlyOneTranscript)) {
        $onlyOneTranscript = 0;
    }

    while(defined($line = <$file>)) {
	@entry = ();
	my $exon;
	@entry = split(/\s+/,$line);
#	print "%%%% ", $line;
	if($line =~ /^\#/ || $line =~ /^\s*$/) {
#	    if($line =~ /^\#.+\>(\S+)\s+/) {
#		$contigKey = $1;
#	    }
	    if($#entry > 0 && $entry[1] =~ /SEQ\:/) { # hmmgene
		$counter++;
	    }
	    next;
	}
	# GFF goes like this : <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments], [attributes] usually have gene "gene_id" as first attribute
	if($entry[2] =~ /gene/i) {
	    $counter++;
	}
	if(defined($entry[8])) {
	    if($entry[8] =~ /gene_id/ || $entry[8] =~ /Sequence/) { #twinscan and others
		$entry[9] =~ /^\"?([^\"]+)\"?\;?/;
		$geneId = ($uniq_id ? $1 : $entry[0].".".$1);
		if(defined($entry[10]) && $entry[10] =~ /transcript_id/) {
		    $entry[11] =~ /^\"?([^\"]+)\"?\;?/;
		    $transcriptId = ($uniq_id ? $1 : $entry[0].".".$1);
		}
		else {
		    $transcriptId = $geneId;
		}
	    }
	    elsif($entry[8] =~ /transcript_id/) { #augustus, encode
		$entry[9] =~ /^\"?([^\"]+)\"?\;?/;
		$transcriptId = $1;
		if(defined($entry[10]) && $entry[10] =~ /gene_id/) {
		    $entry[11] =~ /^\"?([^\"]+)\"?\;?/;
		    $geneId = $1;
		}
		else {
		    $geneId = $transcriptId;
		}

		$transcriptId =~ /^\"?([^\"]+)\"?\;?/;
		$transcriptId = ($uniq_id ? $1 : $entry[0].".".$1);
		$geneId =~  /^\"?([^\"]+)\"?\;?/;
		$geneId = ($uniq_id ? $1 : $entry[0].".".$1);
	    }
	    elsif($entry[8] =~ /trans.+\=[^\|]*\|?([^\";]+)/) { #genezilla 
		$geneId = ($uniq_id ? $1 : $entry[0].".".$1);
		$transcriptId = ($uniq_id ? $1 : $entry[0].".".$1);
	    }
	    elsif($entry[8] =~ /(\S+)parse\:cds_(\d+)/i) {
		$geneId = ($uniq_id ? $counter.$1 : $entry[0].".".$counter.$1);
		$transcriptId = ($uniq_id ? $counter.$1 : $entry[0].".".$counter.$1);
	    }
	    elsif($entry[8] =~ /^ID\=([^\;]+)/i) {
		$geneId = $1;
		$transcriptId = $1;
	    }
	    elsif($entry[8] =~ /$entry[0]/) {
		$entry[8] =~  /^\"?([^\"]+)\"?\;?/;
		$geneId = $1;
		$transcriptId = $1;
	    }
	    else { #all others
		$entry[9] =~ /^\"?([^\"]+)\"?\;?/;
		$geneId = $1;
		$transcriptId = $1;
		
	    }
	}
	else {
	    $entry[0] =~ /^\"?([^\"]+)\"?\;?/;
	    $geneId = $1;
	    $transcriptId = $1;
	}
	$entry[0] =~ /^(\S+)(\.gb\.fa)?/;
	$contigKey = $1;
	
	if($oldGeneId ne "" && $geneId ne $oldGeneId) { # We have changed gene
	    # insert all the alternative transcripts
#	    print STDERR "**** $line\"$transcriptId\" \"$oldGeneId\" \"$geneId\" $contigKey $entry[2]\n"; 
	    GenUtils::insertGFFGene($genes{$source}, \%exons, \%utrExons, \%lengths, \%strand, \%noStart, \%noStop, \%phase, $onlyOneTranscript, $contig, $oldGeneId, $type, $addStart, $preserve_phases);
	    %exons = ();
	    %utrExons = ();
	    %noStart = ();
	    %lengths = ();
	    %noStop = ();
	    %strand  = ();
	    %phase = ();
	}

	$oldGeneId = $geneId;
	
	if(!defined($exons{$transcriptId})) {
	    @{$exons{$transcriptId}} = ();
	    @{$utrExons{$transcriptId}} = ();
	    $noStop{$transcriptId} = 1;
	    $noStart{$transcriptId} = 1;
	    $phase{$transcriptId} = 0;
	}

	if($entry[2] =~ /^single/i) {
	    $noStart{$transcriptId} = 0; $noStop{$transcriptId} = 0;
	}
	elsif($entry[2] eq "start_codon" || $entry[2] =~ /^init/i || $entry[2] =~ /first/i) {
	    $noStart{$transcriptId} = 0;
	    if($entry[2] eq "start_codon") {
		next;
	    }
	    $addStart = 0;
	}
	elsif($entry[2] eq "stop_codon" || $entry[2] =~ /^final/i || $entry[2] =~ /terminal/i || $entry[2] =~ /last/i) {
	    $noStop{$transcriptId} = 0;
	    if($entry[2] eq "stop_codon") {
		next;
	    }
	}
	elsif($entry[2] !~ /est_match/i && $entry[2] !~ /cds/i && $entry[2] !~ /exon/i && $entry[2] !~ /internal/i && $entry[2] !~ /five_prime_UTR/i && $entry[2] !~ /three_prime_UTR/i && $entry[2] !~ /UTR/i ) {
	    next;
	}
	
	if($entry[7] eq '.' && $entry[2] !~ /UTR/i) {
	    print STDERR "Warning: $geneId|$transcriptId contains a nonUTR segment with no phase!\n";
	}

	$numExons++;
	
	my $begExon = $entry[3];
	my $endExon = $entry[4];
# check if coordinates are negative. skip if that's the case
	if($begExon < 1 || $endExon < 1) {
	    next;
	}
	
	$contig = undef;
	if(defined($contigs)) {
	    $contig = $contigs->{$contigKey};
	}
	if(defined($contig)) {
          if($begExon > $contig->{'Len'} || $endExon > $contig->{'Len'}) {
              next;
          }                
	}
	
	$strand{$transcriptId} = ($entry[6] eq "-" ? STRAND_COMP:STRAND_FWD);
	$source = $entry[1];	
	if(!defined($genes{$source})) {
	    %{$genes{$source}} = ();
	}
	if($strand{$transcriptId} == STRAND_COMP) {
	    my $tmp = $begExon;
	    $begExon = $endExon;
	    $endExon = $tmp;
	}
	
	$exon =  [$contigKey, $begExon, $endExon];
	if($entry[2] =~ /five_prime_UTR/i || $entry[2] =~ /three_prime_UTR/i || $entry[2] =~ /UTR/i) {
	    $exons = $utrExons{$transcriptId};
	}
	else {
	    $exons = $exons{$transcriptId};
	}
#	print "exon ", join(" ", @$exon), " " , $strand{$transcriptId}, "\n";
# check if exon is already inserted
	my $alreadyInserted = 0;
	for(my $i = 0; $i < scalar(@{$exons}); $i++) {
	    if($exon->[1] == $exons->[$i]->[1] || $exon->[2] == $exons->[$i]->[2]) {
		if(abs($exon->[2] - $exon->[1]) < abs($exons->[$i]->[2] - $exons->[$i]->[1])) { # it happens when $exons[$i] is not totally coding
		    $exons->[$i]->[1] = $exon->[1];
		    $exons->[$i]->[2] = $exon->[2];
		}
		$alreadyInserted = 1;
		last;
	    }
	}
	if($alreadyInserted) {
          next;
	}
	my $j = 0;

	if($strand{$transcriptId} == STRAND_COMP) {	    
	    for($j = 0; $j < scalar(@{$exons}); $j++) {
		if($exon->[1] >= $exons->[$j]->[1]) {
		    last;
		}
	    }

	    @{$exons} = (@{$exons}[0 .. $j - 1], $exon, 
			 @{$exons}[$j .. $#$exons]);
	}
	else {
	    for($j = 0; $j < scalar(@{$exons}); $j++) {
		if($exon->[1] <= $exons->[$j]->[1]) {
		    last;
		}
	    }

	    @{$exons} = (@{$exons}[0 .. $j - 1], $exon,
			 @{$exons}[$j .. $#$exons]);
	}
	
	if($entry[2] =~ /three_primer_UTR/i || $entry[2] =~ /five_primer_UTR/i || $entry[2] =~ /UTR/i) {
	    next;
	}
	
	if($strand{$transcriptId} == STRAND_COMP) {
	    if(scalar(@{$exons}) == 1 || $exon->[1] > $exons->[-1]->[1]) {
		$phase{$transcriptId} = ($entry[7] eq '.' ? 0 : (3 - $entry[7]) % 3);
	    }
	}
	else {
	    if(scalar(@{$exons}) == 1 || $exon->[1] < $exons->[0]->[1]) {
		$phase{$transcriptId} = ($entry[7] eq '.' ? 0 : (3 - $entry[7]) % 3);
	    }
	}

	$lengths{$transcriptId} += abs($exon->[2] - $exon->[1]);
	
    }
    
    # Inserting last gene
    GenUtils::insertGFFGene($genes{$source}, \%exons, \%utrExons, \%lengths, \%strand, \%noStart, \%noStop, \%phase, $onlyOneTranscript, $contig, $oldGeneId, $type, $addStart, $preserve_phases);

    print STDERR "Total number of exons processed = ", $numExons, "\n";
    return \%genes;
}

sub readLocsFormat2Array {
    my ($file, $includeSeq, $type, $ref_line) = @_;
    my @genes = ();
    my $id;
    my $gid;
    my $tid;
    $$ref_line = <$file>;
    $includeSeq = 0 unless defined($includeSeq);

    while(defined($$ref_line) && ($$ref_line =~ /^>(\S+)\s+(.+)\s+?(\d+)?b?p?$/)) {
        my $locs = ();
        my @locs = ();
        my $seq = "";
        my @seq = ();
        $id = $1;
        my $location = $2;
        my $len = $3;

        if($id =~ /(\S+)\|([^\|]+)/) {
            $gid = $1;
            $tid = $2;
        }
        else {
            undef $gid;
            undef $tid;
        }
        
        my ($noStart, $phase, $fiveUtrExons, $exons, $threeUtrExons, $noStop) = @{&readContigLocation($location)};
#	print "detail $id \"$location\" $noStart, $phase, $#$fiveUtrExons, $#$exons, $#$threeUtrExons, $noStop\n";
        if(defined($includeSeq) && $includeSeq) {
            while(defined($$ref_line = <$file>) && $$ref_line !~ /^>/) {
                $$ref_line =~ s/\s+//g;
                push(@seq, $$ref_line);
            }
            $seq = join("", @seq);
        }
        else {
            $$ref_line = <$file>;
        }

        my $thisExons = $exons;

        if(!scalar(@$exons)) {
            $thisExons = scalar(@$fiveUtrExons) ? 
                $fiveUtrExons :
                $threeUtrExons;
        }
        
#	    print "\"$id $gid $tid\"\n";
        my $gene = Gene->new('Id' => $id, 'tId' => $tid, 'gId' => $gid, '5UTRExons' => $fiveUtrExons, 'Exons' => $exons, '3UTRExons' => $threeUtrExons, 'Contig' => undef, 'Seq' => $seq, 'Strand' => ($thisExons->[-1]->[2] > $thisExons->[0]->[1]) ? STRAND_FWD : STRAND_COMP, 'NoStart' => $noStart, 'NoStop' => $noStop, 'Phase' => $phase);
#	print STDERR $gene->{'Id'}, " ", scalar(@{$gene->{'Exons'}}),"\n";
	push(@genes, [$id, $gene]);
    }

    return \@genes;

}

sub readLocsFormat {
    my ($file, $includeSeq, $type, $contigs) = @_;
    my %genes = ();
    my $line = <$file>;
    $includeSeq = 0 unless defined($includeSeq);

    while(1) {
	my $gene = readLocsFromString($line, $contigs);
	if(!defined($gene)) {
	    last;
	}

        if(defined($includeSeq) && $includeSeq) {
	    my @seq = ();
            while(defined($line = <$file>) && $line !~ /^>/) {
                $line =~ s/\s+//g;
                push(@seq, $line);
            }
	    $gene->{'Seq'} = join("", @seq);
        }
        else {
            $line = <$file>;
        }

	my $id = $gene->{'Id'};
	my $contig = $gene->{'Contig'};
        if(defined($genes{$id})) {
            warn "Repeated entry $id\n";
        }

        $genes{$id} = $gene;
        if(defined($contig)) {
            $contig->{'Features'}->{$type}->{$id} = $gene;
        }

    }

    return \%genes;
}

sub readLocsFromString {
    my ($line, $contigs) = @_;
    my $id;
    my $gid;
    my $tid;

    if(defined($line) && ($line =~ /^>(\S+)\s+(.+)\s+?(\d+)?b?p?$/)) {
        my $locs = ();
        my @locs = ();
        $id = $1;
        my $location = $2;
        my $len = $3;

        if($id =~ /(\S+)\|(\[^\|]+)/) {
            $gid = $1;
            $tid = $2;
        }
        else {
            undef $gid;
            undef $tid;
        }
        
        my ($noStart, $phase, $fiveUtrExons, $exons, $threeUtrExons, $noStop) = @{&readContigLocation($location)};
        my $contig = undef;
        my $thisExons = $exons;

        if(!scalar(@$exons)) {
            $thisExons = scalar(@$fiveUtrExons) ? 
                $fiveUtrExons :
                $threeUtrExons;
        }
        
        if(defined($contigs)) {

            defined($contigs->{$thisExons->[0]->[0]}) or die "Cannot find contig \"$thisExons->[0]->[0]\" in $id\n";
            $contig = $contigs->{$thisExons->[0]->[0]};

        }
#	print STDERR "\"$id $tid\"\n";
#	print STDERR "$id \"$location\" $noStart, $phase.. scalars4exons ", scalar(@$fiveUtrExons), " ", scalar(@$exons), " ", scalar(@$threeUtrExons), " $noStop\n";

        return Gene->new('Id' => $id, 'tId' => $tid, 'gId' => $gid, '5UTRExons' => $fiveUtrExons, 'Exons' => $exons, '3UTRExons' => $threeUtrExons, 'Contig' => $contig, 'Seq' => "", 'Strand' => ($thisExons->[-1]->[2] > $thisExons->[0]->[1]) ? STRAND_FWD : STRAND_COMP, 'NoStart' => $noStart, 'NoStop' => $noStop, 'Phase' => $phase);
    }
    else { 
	return undef;
    }
}


sub readWeightedFastaFormat {
    my ($file, $includeSeq, $type, $contigs) = @_;
    my %genes = ();
    my $id;
    my $gid;
    my $tid;
    my $line = <$file>;
    $includeSeq = 0 unless defined($includeSeq);
    
    while(defined($line) && ($line =~ /^>(\S+)\s+([^\s]+)\s+(.+)$/)) {
#	print STDERR $line;
        my $locs = ();
        my @locs = ();
        my $seq = "";
        my @seq = ();
        $id = $1;
        my $location = $2;
	my @weight = split(/\;/, $3);

        if($id =~ /(\S+)\|([^\|]+)/) {
            $gid = $1;
            $tid = $2;
        }
        else {
            undef $gid;
            undef $tid;
        }
        
        my ($noStart, $phase, $fiveUtrExons, $exons, $threeUtrExons, $noStop) = @{&readContigLocation($location)};

	for(my $i = 0; $i < scalar(@$exons); $i++) {
	    push(@{$exons->[$i]}, $weight[$i]);
	}

        if(defined($includeSeq) && $includeSeq) {
            while(defined($line = <$file>) && $line !~ /^>/) {
                $line =~ s/\s+//g;
                push(@seq, $line);
            }
            $seq = join("", @seq);
        }
        else {
            $line = <$file>;
        }

        my $contig = undef;
        my $thisExons = $exons;

        if(!scalar(@$exons)) {
            $thisExons = scalar(@$fiveUtrExons) ? 
                $fiveUtrExons :
                $threeUtrExons;
        }
        
        if(defined($contigs)) {

            defined($contigs->{$thisExons->[0]->[0]}) or die "Cannot find contig \"$thisExons->[0]->[0]\" in $id\n";
            $contig = $contigs->{$thisExons->[0]->[0]};

        }

        my $gene = Gene->new('Id' => $id, 'tId' => $tid, 'gId' => $gid, '5UTRExons' => $fiveUtrExons, 'Exons' => $exons, '3UTRExons' => $threeUtrExons, 'Contig' => $contig, 'Seq' => $seq, 'Strand' => ($thisExons->[-1]->[2] > $thisExons->[0]->[1]) ? STRAND_FWD : STRAND_COMP, 'NoStart' => $noStart, 'NoStop' => $noStop, 'Phase' => $phase);
        if(defined($genes{$id})) {
            warn "Repeated entry $id\n";
        }
#	print STDERR $gene->{'Id'}, " ", scalar(@{$gene->{'Exons'}}),"\n";

        $genes{$id} = $gene;
        if(defined($contig)) {
            $contig->{'Features'}->{$type}->{$id} = $gene;
        }
    }

    return \%genes;
}


### Read exons in USCS format, offset is substracted for each coordinate
sub readUCSCFormat {
    my ($file, $includeSeq, $offset) = @_;
    my $line = "";
    my ($contig, $id, $exon, $oldId, $begin, $end, $revcomp) = ("","","","",0,0,0);
    my @exons = ();
    my %res = ();
    my $seq = "";

    if(!defined($offset)) {
	$offset = 0;
    }
    $line = <$file>;
    while(defined($line) && $line =~ />(\S+)_(\d+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+(\(reverse complemented\))?/) {
        my @exonSeq = ();
        $id = $1;
        $exon = $2;
        $contig = $3;
        $begin = $4 - $offset;
        $end = $5 - $offset;
	$revcomp = 0;

        if(defined($6) && $6 =~ /reverse/) {
            $revcomp = 1;
        }

        if(($oldId ne "" && $oldId ne $id) || !defined($line)) { ### change id
            my $thisSeq = $seq;
            my @thisExon = ();
            my $loc;
	    foreach $loc (@exons) {
                push(@thisExon, $loc);
            }
#	    print "$oldId $thisSeq\n";
	    my $gene = Gene->new('Id' => $oldId, 'Seq' => $thisSeq, 'Exons' => \@thisExon, 'Strand' => ($loc->[2] >= $loc->[1] ? STRAND_FWD : STRAND_COMP));
            $res{$oldId} = $gene;
            $seq = "";
            @exons = ();
        }

        if(defined($includeSeq) && $includeSeq) {
            while(defined($line = <$file>) && $line !~ />(\S+)_(\d+)\s+\S+\s+\S+(\d+)\s+(\d+)\s+(\(reverse complemented\))?/) {
                $line =~ s/\s+//g;
                push(@exonSeq, $line);
            }
        }
        else {
            $line = <$file>;
        }

        if($revcomp) {
            my $tmp = $end;
            $end = $begin;
	    $begin = $tmp;
            $seq = join("",@exonSeq).$seq;
            unshift (@exons, [$contig, $begin, $end]); 
        }
        else {
            $seq = $seq.join("",@exonSeq);
            push(@exons, [$contig, $begin, $end]);
        }
        $oldId = $id;
    }
    return \%res;
}

sub readLocExons {
    my($location, $weights) = @_;
    my @exons = ();
    my @locs = ();
    my @weights = ();
    $weights = "" unless defined($weights);

    if($location =~ /\;/) {
	@locs = split(/;/, $location);
    }
    else {
        push(@locs, $location);
    }

    (scalar(@locs) > 0) or die $location."\n";
    
    if($weights =~ /\;/) {
	@weights = split(/\;/, $weights);
    }
    elsif($weights ne "") {
	for(my $i = 0; $i <= $#locs; $i++) {
	    push(@weights, $weights/scalar(@locs));
	}
    }

    for(my $i = 0; $i <= $#locs; $i++) {
	my $exon = $locs[$i];
        ($exon =~ /^\s*(\S+)_(\d+)_(\d+)\s*/) or die "bad exon $location $weights\n";
        if($exon =~ /^\s*(\S+)_(\d+)_(\d+)\s*/) {
	    if(scalar(@weights)) {
		($#weights == $#locs) || die "wrong number of weights for location $location\n";
		push(@exons ,[$1,$2,$3, $weights[$i]]);
	    }
	    else {
		push(@exons ,[$1,$2,$3]);
	    }
        }
    }
    return \@exons;
}


sub readContigLocation {
    my ($location) = @_;
    my @locs = ();
    my $fiveUtrExons;
    my $exons;
    my $threeUtrExons;
    my $exon;
    my ($noStart, $phase, $noStop) = (0, 0, 0);
    my $weights = "";

    if($location =~ /^([^\s]+)\s+([\d,\],\[,\.,\,\-;]+)$/) { # weighted locations
	$location = $1;
	$weights = $2;
    }

    if($location =~ /([^\[]+)\[(\S*)\s*$/) {
        $location = $2;
	my $fivePlocation = $1;
	my $fivePweights = "";
	if($weights =~ /([^\[]+)\[(\S*)\s*$/) {	
	    $fivePweights = $1;
	    $weights = $2;
	}
        $fiveUtrExons = &readLocExons($fivePlocation, $fivePweights);
    }
    else {
        @{$fiveUtrExons} = ();
    }

    if($location =~ /([^\]]*)\](\S+)\s*$/) {
        $location = $1;
	my $threePlocation = $2;
	my $threePweights = "";
	if($weights =~ /([^\]]*)\](\S+)\s*$/) {
	    $weights = $1;
	    $threePweights = $2;
	}	

        $threeUtrExons = &readLocExons($threePlocation, $threePweights);
    }
    else {
        @{$threeUtrExons} = ();
    }

    if($location =~ /(\d+)?\<(\S*)\s*$/) {
        if(defined($1)) {
            $phase = $1;
        }
        $location = $2;
	$noStart = 1; # partial first exon
    }
    
    if($location =~ /^(\S*)\>/) {
        $location = $1;
        $noStop = 1; # partial last exon 
    }

    if($location =~ /^([^\s]+)$/) {	
        $exons = &readLocExons($location, $weights);
    }
    else {
        @{$exons} = ();
    }

    return [$noStart, $phase, $fiveUtrExons, $exons, $threeUtrExons, $noStop];
}

#### OBSOLETE FUNCTION #####
sub findPossibleGeneSuccesors {
    my($genes, $indexes, $maxDistance) = @_;
    my $i;
    my @indGenes = ();
    my @results = ();
    for($i = 0; $i < scalar(@$indexes); $i++) {
	push(@indGenes, $genes->[$indexes->[$i]]);
    }
    my @sortedIndGenes = sort Gene::compareTo (@indGenes);
    my ($lowerBound, $upperBound) = (0,0);
    if($sortedIndGenes[0]->{'Strand'} == STRAND_FWD) {
	$lowerBound = $sortedIndGenes[0]->{'Exons'}->[-1]->[2];
    }
    else {
	$lowerBound = $sortedIndGenes[0]->{'Exons'}->[0]->[1];
    }
    if($sortedIndGenes[-1]->{'Strand'} == STRAND_FWD) {    
	$upperBound = $sortedIndGenes[-1]->{'Exons'}->[-1]->[2] + $maxDistance;
    }
    else {
	$upperBound = $sortedIndGenes[-1]->{'Exons'}->[0]->[1] + $maxDistance;
    }
    $i = 0;
    while(1) {
	if($genes->[$i]->{'Strand'} == STRAND_FWD) {
	    if($genes->[$i]->{'Exons'}->[0]->[1] > $lowerBound) {
		if($genes->[$i]->{'Exons'}->[0]->[1] < $upperBound) {
		    push(@results, $i);
		}
		else { 
		    last;
		}
	    }
	}
	else {
	    if($genes->[$i]->{'Exons'}->[-1]->[2] > $lowerBound) {
		if($genes->[$i]->{'Exons'}->[-1]->[2] < $upperBound) {	
		    push(@results, $i);		
		}    
		else {
		    last;
		}
	    }
	}
	$i++;
    }
    return \@results;
}

#Writes locs and associated contigs in Glimmer format which we have copied and
# pasted from command line
#    <exon_file>  is a file with the exon coordinates relative to the sequences
#    contained in the <mfasta_file>; different genes are separated
#    by a blank line; I am assuming a format like below:
#
#                 seq1 5 15
#                 seq1 20 34
#
#                 seq1 50 48
#                 seq1 45 36

#                 seq2 17 20

#                 In this example seq1 has two genes: one on the direct strand 
#                 and another one on the complementary strand

#                 The partial genes can be specified as in the following example:
#    
#                 seq2 <100 >234
#
#                 seq3 <1  100
#                 seq3 102 >105

sub writeGlimmerFormat {
    my($locs, $exons_file) = @_;
    my @keysList = keys %$locs;
    my ($i, $exon);

    for(my $count = 0; $count < scalar(@keysList); $count++)  {
	my $gene = $locs->{$keysList[$count]};
	for($i = 0; $i < scalar(@{$gene->{'Exons'}}) ; $i++) {

	    $exon = $gene->{'Exons'}->[$i];
	    print $exons_file $exon->[0], " ";
	    if($i  == 0 && $gene->{'NoStart'}) {
		print $exons_file "<";
	    }
	    print $exons_file $exon->[1], " ";
	    if($i == scalar(@ { $gene->{'Exons'}}) - 1 && $gene->{'NoStop'}) {
		print $exons_file ">";
	    }
	    print $exons_file $exon->[2], "\n";
	}
	print $exons_file "\n";
    }
}

#Format is something like this
#GlimmerM (Version 3.0)
#Sequence name: D49493
#Sequence length: 17286 bp
#
#Predicted genes/exons
#
#Gene Exon Strand  Exon            Exon Range      Exon
#   #    #         Type                           Length
#
#   1    1  +  Initial         525        557       33
#   1    2  +  Internal       1625       1740      116
#   1    3  +  Internal       3227       4235     1009
#   1    4  +  Internal       4332       4338        7
#   1    5  +  Internal       4392       4502      111

sub readGlimmerFormat {
    my($file) = @_;
    my $i  = 0;
    my ($line, $id);
    my %genes = ();
    while(defined($line = <$file>)) {
	if(defined($line) && $line !~ /Sequence name: (\S+)$/) {
	    $id = $1;
	}
	$i = 0;
	while(defined($line) && (($line =~ /Initial/ && $line =~ /\+/) || ($line =~ /Terminal/ && $line =~ /\-/) || $line =~ /Single/)) {
	    my @exons = ();
	    my $strand;
	    while(1) {
		my @line = split(/\s+/, $line);
		$strand = ($line[3] eq "+") ? STRAND_FWD : STRAND_COMP;
		if($strand == STRAND_COMP) {
		    unshift(@exons, [$id, $line[6], $line[5]]);
		}
		else {
		    push(@exons, [$id, $line[5], $line[6]]);
		}
		if(($line =~  /Terminal/ &&  $line =~ /\+/) || ($line =~ /Initial/ && $line =~ /\-/) || $line =~ /Single/) {
		    last;
		}
		$line = <$file>;
	    } 
#	    print $id.$i," ",join(" ", @exons), "\n";
	    $line = <$file>; $line = <$file>;
	    my $gene = Gene->new('Id' => "$id.".$i, 'Exons' => \@exons, 'Seq' => "", 'Strand' => $strand);
	    $genes{"$id.".$i} = $gene;
	    $i++;
	}	
    } 

    return \%genes;
}

####################################################################################
# sub readGenScanFormat
# exon coordinates follow GTF format. Here is an example
#  Gn.Ex Type S .Begin ...End .Len Fr Ph I/Ac Do/T CodRg P.... Tscr..
# ----- ---- - ------ ------ ---- -- -- ---- ---- ----- ----- ------
#
# 1.01 Intr +   3620   3731  112  0  1   58   44    87 0.271   0.74
# 1.02 Intr +   5434   5516   83  1  2  113   82   103 0.986  11.48
# 1.03 Intr +   5810   5911  102  0  0  101  107   141 0.998  17.45
# 1.04 Intr +   6665   6739   75  0  0  107   94    98 0.998  11.79
# 1.05 Term +   7637   8415  779  0  2   76   42   816 0.997  69.33
# 1.06 PlyA +   8545   8550    6                               1.05
#
# 2.02 PlyA -   8604   8599    6                               1.05
# 2.01 Term -   9243   9166   78  0  0   78   45    96 0.326   2.06
#
####################################################################################

sub readGenScanFormat {
    my($file, $type, $contigs, $preserve_phases) = @_;
    
    $preserve_phases = 0 if !defined($preserve_phases);

    my %genes = ();
    my $numExons = 0;
    my $strand = STRAND_FWD;
    my @line = ();
    my $line = <$file>;
    my $id;
    my @exons = ();
    my $geneNum;
    my $phase;
    my $hasStart;
    my $hasStop;
    my $contig;
    while(defined($line)) {
        while(defined($line) && $line !~ /Sequence\s+(\S+)\s+\:/) { $line = <$file>;}
        defined($line) || last;
        $line =~ /Sequence\s+(\S+)\s+\:/;
        $id = $1;
        while(defined($line) && $line !~ /Gn\.Ex/) { $line = <$file>;}
        @exons = ();
        $geneNum = 1;
        $phase = -1;
        $hasStart = 0;
        $hasStop = 0;
        my $gene;

        while(defined($line = <$file>) && $line !~ /^Sequence/) {
            my $cdsRegion = 0;
            $line =~ /^\s*(\S.+\S)\s*$/;
            $line = $1;
            @line = split(/\s+/, $line);
#	    print  scalar(@line),  " ", $line[1]," ", $line;
            (scalar(@line) == 13 && ($line[2] eq "+" || $line[2] eq "-")) || next;
            
            if($line[1] eq "Intr" || $line[1] eq "Term" || $line[1] eq "Init" || $line[1] eq "Sngl") {
                $cdsRegion = 1;
            }
            ($line[0] =~ /(\d+)\.(\d+)/) || next;
            
            if($1 != $geneNum && scalar(@exons) > 0) { # a new gene
                my @tmp = @exons;
#		print "----> $id $geneNum $strand ", scalar(@tmp),"\n";
                $contig = undef;
                $numExons += scalar(@tmp);
                if(defined($contigs)) {
                    $contig = $contigs->{$exons[0]->[0]};
                }
                $gene = Gene->new('Id' => "$id.".$geneNum, 'Exons' => \@tmp, 'Seq' => "", 'Strand' => $strand, 'NoStart' => !$hasStart, 'NoStop' => !$hasStop, 'Phase' => $phase, 'Contig' => $contig);
		if(defined($contig)) {
		    $contig->{'Features'}->{$type}->{"$id.".$geneNum} = $gene;
		    $gene->checkStop();
		    $gene->computePhase() if !$preserve_phases;
		    $gene->checkStart();
		}
                $genes{"$id.".$geneNum} = $gene;
                $geneNum = $1;
                $hasStart = 0;
                $hasStop = 0;
                $phase = -1;
                @exons = ();
            }
            
            if($line[1] eq "Init" || $line[1] eq "Sngl") {
                $hasStart = 1;
                $phase = 0;
            }
            if($line[1] eq "Intr" && $phase == -1) {
                $phase = $line[7];
            }
            if ($line[1] eq "Term" || $line[1] eq "Sngl") {
                $hasStop = 1; 
            }
            $strand = ($line[2] eq "+") ? STRAND_FWD : STRAND_COMP;
            
            if($cdsRegion) {
                if($strand == STRAND_COMP) {
                    unshift(@exons, [$id, $line[3], $line[4]]);
                }
                else {
                    push(@exons, [$id, $line[3], $line[4]]);
                }
            }
        }

        if(scalar(@exons) > 0) {
#	    print "<$id $geneNum $strand ", scalar(@exons),"\n";
            $contig = undef;
            if(defined($contigs)) {
                $contig = $contigs->{$exons[0]->[0]};
            }
            $numExons += scalar(@exons);
            $gene = Gene->new('Id' => "$id.".$geneNum, 'Exons' => \@exons, 'Seq' => "", 'Strand' => $strand, 'NoStart' => !$hasStart, 'NoStop' => !$hasStop, 'Phase' => $phase, 'Contig' => $contig);
	    if(defined($contig)) {
		$contig->{'Features'}->{$type}->{"$id.".$geneNum} = $gene;
		$gene->checkStop();
		$gene->computePhase() if !$preserve_phases;
		$gene->checkStart();
	    }
            $genes{"$id.".$geneNum} = $gene;
        }

    }
    print STDERR "Total number of exons processed = ", $numExons, "\n";    
    return \%genes;
}

###############################################################################
# sub writeGTFFormat 
###############################################################################

sub writeGTFFormat {
    my($locs, $gff_file, $contigs) = @_;
    my @keysList = keys %$locs;
    for(my $count = 0; $count < scalar(@keysList); $count++)  {
	my $gene = $locs->{$keysList[$count]};
	writeGeneGTFFormat($gene, $gff_file);
    }
}

sub writeGeneGTFFormat {
    my($gene, $gff_file, $tr_field) = @_;
    my $i;
    my $exon;
    my $frame = 0;
    my $geneLength = $gene->{'Phase'};
    my $strand = "+";
    my $score;
    
    my $geneId = $gene->{'Id'};
    my $transcriptId = $gene->{'Id'}.".".$gene->{'Strand'};
    my $tr_field_printed = 0;
    
    if($gene->{'Id'} =~ /(\S+)\|([^\|]+)/) {
	$geneId = $1;
	$transcriptId = $2;
    }

    for($i = 0; $i < scalar(@{$gene->{'5UTRExons'}}) ; $i++) {
	$exon = $gene->{'5UTRExons'}->[$i];
	my $begExon = $exon->[1];
	my $endExon = $exon->[2];
	if($gene->{'Strand'} == STRAND_COMP) {
	    my $tmp = $begExon;
	    $begExon = $endExon;
	    $endExon = $tmp;
	    $strand = "-";
	}
	$score = defined($exon->[3]) ? "$exon->[3]" : ".";
	print $gff_file "$exon->[0]\tCRAIG\t5UTR\t$begExon\t$endExon\t$score\t$strand\t.\tgene_id \"$geneId\"; transcript_id \"$transcriptId\";";
	if(defined ($tr_field)) {
	    print $gff_file "\t", $tr_field;
	    $tr_field_printed = 1;
	}
	print "\n";
    }
    
    $strand = "+";
    for($i = 0; $i < scalar(@{$gene->{'Exons'}}) ; $i++) {
	$exon = $gene->{'Exons'}->[$i];
	$frame = (3 - ($geneLength % 3)) % 3;
	$score = defined($exon->[3]) ? "$exon->[3]" : ".";
	my $begExon = $exon->[1];
	my $endExon = $exon->[2];
	my $startBeg = $begExon;
	my $startEnd = $begExon + 2;
	my $stopBeg = $endExon + 1;
	my $stopEnd = $endExon + 3;
	if($gene->{'Strand'} == STRAND_COMP) {
	    my $tmp = $begExon;
	    $startBeg = $begExon - 2;
	    $startEnd = $begExon;
	    $stopBeg = $endExon - 3;
	    $stopEnd = $endExon - 1;
	    $begExon = $endExon;
	    $endExon = $tmp;
	    $strand = "-";
	}
	
	if($i == 0 && !$gene->{'NoStart'}) {
	    print $gff_file "$exon->[0]\tCRAIG\tstart_codon\t$startBeg\t",$startEnd, "\t$score\t$strand\t0\tgene_id \"$geneId\"; transcript_id \"$transcriptId\";";
	    if(!$tr_field_printed && defined ($tr_field)) {
		print $gff_file "\t", $tr_field;
		$tr_field_printed = 1;
	    }
	    print "\n";
	}
	
	print $gff_file "$exon->[0]\tCRAIG\tCDS\t$begExon\t$endExon\t$score\t$strand\t$frame\tgene_id \"$geneId\"; transcript_id \"$transcriptId\";";
	if(!$tr_field_printed && defined ($tr_field)) {
	    print $gff_file "\t", $tr_field;
	    $tr_field_printed = 1;
	}
	print "\n";
	
	if($i == scalar(@ { $gene->{'Exons'}}) - 1 && !$gene->{'NoStop'}) {
	    print $gff_file "$exon->[0]\tCRAIG\tstop_codon\t", $stopBeg, "\t", $stopEnd, "\t$score\t$strand\t0\tgene_id \"$geneId\"; transcript_id \"$transcriptId\";\n";
	}
	$geneLength += abs($exon->[2] - $exon->[1]) + 1;
    }
    
    $strand = "+";
    for($i = 0; $i < scalar(@{$gene->{'3UTRExons'}}) ; $i++) {
	$exon = $gene->{'3UTRExons'}->[$i];
	my $begExon = $exon->[1];
	my $endExon = $exon->[2];
	if($gene->{'Strand'} == STRAND_COMP) {
	    my $tmp = $begExon;
	    $begExon = $endExon;
	    $endExon = $tmp;
	    $strand = "-";
	}
	$score = defined($exon->[3]) ? $exon->[3] : ".";
	print $gff_file "$exon->[0]\tCRAIG\t3UTR\t$begExon\t$endExon\t$score\t$strand\t.\tgene_id \"$geneId\"; transcript_id \"$transcriptId\";";
	if(!$tr_field_printed && defined ($tr_field)) {
	    print $gff_file "\t", $tr_field;
	    $tr_field_printed = 1;
	}
	print "\n";
    }
}


sub displayContigs {
    my($list, $file) = @_;
    my($id);
    foreach $id (keys %$list) {
	print $file ">$list->{$id}->{'Id'}\n$list->{$id}->{'Seq'}\n";
    }
}

sub displayGeneList {
    my($list, $file, $contigs, $onlyHeader, $beg, $end) = @_;
    my($id);
    my $count = 0;
    my @keysList = keys %$list;
    $beg = 0 if (!defined($beg) || $beg < 0);
    $end = scalar(@keysList) if (!defined($end) || $end > scalar(@keysList));
    for($count = $beg; $count < $end; $count++)  {
        $id = $keysList[$count];
        $list->{$id}->dump($file, $contigs, $onlyHeader);
    }
}

sub readGenBankLocation {
    my ($location, $cid, $offset, $length) = @_;
    my @exons = ();
    my @tmp_exons = ();
    my ($noStart, $noStop) = (0, 0);
    my ($i, $beg, $end);
    my $revStrand = 0;
    my $nostart_str;
    my $nostop_str;

    if($location =~ /complement\((.+)\)$/) {
	$location = $1;
	$revStrand = 1;
    }
#	print $line;
    
    if($location =~ /join\((.+)\)$/) {
	@tmp_exons = split(/\,/,$1);
    }
    else {
	push(@tmp_exons, $location);
    }
    
    for($i = 0; $i < scalar(@tmp_exons); $i++) {
	my $exon = $tmp_exons[$i];
	$exon =~ s/\s+//g;
	if($exon =~ /^(\<?)(\d+)(\>?)$/) {
	    $nostart_str = $1;
	    $beg = $end = $2;
	    $nostop_str = $3;
	}
	elsif($exon =~ /^(\<?)(\d+)\.\.(\>?)(\d+)$/) {
	    $nostart_str = $1;
	    $beg = $2;
	    $nostop_str = $3;
	    $end = $4;
	}
	else {
	    @exons = ();
	    last;
	}

	if($i == 0 && defined($nostart_str) && $nostart_str eq "<") {
	    $noStart = 1;
	}
	($end >= $beg) || die "$cid $beg $end\n"; # true only for CDNAs
	if($revStrand) {
	    my $tmp = $beg;
	    $beg = $end;
	    $end = $tmp;
	}
	if($i == scalar(@tmp_exons) && defined($nostop_str) && $nostop_str eq ">") {
	    $noStop = 1;
	}

	push(@exons, [$cid, $beg, $end]);
#		print "->", "$cid, $beg, $end\n";
    }
    
    @exons = reverse(@exons) if($revStrand);
    return [$noStart, \@exons, $noStop, $revStrand];
}

#-complete option will only print all the CDNA sequence... not only the coding part. This subroutine reads gene structures, a CDS or mRNA feature is required
sub readGenBankFormat {
    my($complete, $file, $type, $cds_tag, $tr_tag, $contigs, $add_contigs) = @_;
    my($line, $id, $id_add, $gseq, $seq);

    my %genes = ();
    my $length = 0;
    my $cds_only = 0;
    if(!defined($add_contigs)) { $add_contigs = 0;}
    if($tr_tag eq 'null') { 
	$cds_only = 1;
    }

    my $cid = "";
    my %xcripts = ();
    my $feat_space = 0;
    my $tid = "";
    while(defined($line = <$file>)) {
	if($line =~ /^\s*LOCUS\s+(\S+)\s+(\d+)\s*bp/) {
	    %xcripts = ();
	    $cid = $1;
	    if(defined($2)) { 
		$length = $2;
	    }
	}
	elsif($line =~ /^\s*FEATURES\s*/) {
	    # extracting source coordinates
	    $feat_space = 1;
	}
	elsif($line =~ /^\s*ORIGIN\s*/) {
	    $feat_space = 0;
	    # get the source sequence
	    $seq = "";
	    while(defined($line = <$file>) && $line !~ /^\/\//) {
		$line =~ s/\d+//g;
		$line =~ s/\s+//g;
		$seq .= $line;
	    }
		    
	    my $contig = undef;
	    if(defined($contigs)) {
		$contig = $contigs->{$cid};
	    }

	    if(!defined($contig) && $add_contigs) {
		$contig = $contigs->{$cid} = Contig->new('Id' => $cid, 'Seq' => $seq, 'Beg' => 1);
	    }

	    foreach $tid (keys %xcripts) {
		my $cds = $xcripts{$tid}->[0];
		my $gene = $cds;
		
		if(defined($cds)) {
		    if($complete && 
		       ($cds->{'NoStart'} || $cds->{'noStop'})) {
			next;
		    }
		    if(defined($xcripts{$tid}->[1])) {
#			print STDOUT "\n merging $tid \"";
#			$cds->displayHeader(\*STDOUT);
#			print "\"\n and ";
#			$xcripts{$tid}->[1]->displayHeader(\*STDOUT), "\n";
			$cds->addUTRAnnot($xcripts{$tid}->[1]);
		    }				
		}
		elsif(!$cds_only) {
		    $gene = $xcripts{$tid}->[1];
		}

		if(!defined($gene)) {
		    next;
		}

		if(defined($contig)) {
		    $contig->{'Features'}->{$type}->{$gene->{'Id'}} = $gene;
		}
		
		$genes{$gene->{'Id'}} = $gene;
	    }
	    $tid = "";
	}

	if($feat_space == 1) {
	    my $offset = 1;
	    my $roffset = $length;
	    if($line =~ /^\s+source\s+(\d+)\.\.(\d+)/) {
		$offset = $1;
		$roffset = $2;
	    }
	    elsif($line =~ /^\s+gene\s+\S+/) {
		$line = <$file>;
		$line =~ /^\s+\S+\=\"?([^\"]+)\"?/;
#		$gid = $1;
#		print "\ngene $1\n";
	    }
	    elsif(!$cds_only && $line =~ /^\s+$tr_tag\s+(\S+)\s*$/) {
		my $location = $1;
		while(defined($line = <$file>) && $line !~ /^\s+\S+\=/) {
		    chomp $line;
		    $location = $location.$line;
		}
		$location =~ s/\s//g;
		my($noStart, $exons, $noStop, $strand) = @{&readGenBankLocation($location, $cid, $offset, $roffset - $offset + 1)};
		$line =~ /^\s+\S+\=\"?[^\|]*\|?([^\"]+)\"?/;
		my $tid = $1;
		my $mrna = Gene->new('Id' => "$tid", 'Exons' => $exons, 'Contig' => undef, 'Seq' => undef, 'Strand' => $strand, 'NoStart' => $noStart, 'NoStop' => $noStop);
#		print "\nmrna $tid\t$location\n";
#		$mrna->displayHeader(\*STDOUT), "\n";
		if(!defined($xcripts{$tid})) {
		    $xcripts{$tid} = [undef, $mrna];
		}
		else {
		    $xcripts{$tid}->[1] = $mrna;
		}
	    }
	    elsif($line =~ /^\s*$cds_tag\s+(\S+)\s*$/) {
		my $location = $1;
		while(defined($line = <$file>) && $line !~ /^\s+\S+\=/) {
		    chomp $line;
		    $location = $location.$line;
		}
		$location =~ s/\s//g;
		my($noStart, $exons, $noStop, $strand) = @{&readGenBankLocation($location, $cid, $offset, $roffset - $offset + 1)};
		# trying to find encoded by transcript
		$line =~ /^\s+\S+\=\"?[^\|]*\|?([^\"]+)\"?/;
		my $tid = $1;
		while($line !~ /^\s+\/translation\=/) {
		    if($line =~ /encoded by transcript [^\|]*\|?([^\s,\",\;]+)/) {
			$tid = $1;
			last;
		    }
		    $line = <$file>;
		}

		my $cds = Gene->new('Id' => "$tid", 'Exons' => $exons, 'Contig' => undef, 'Seq' => undef, 'Strand' => $strand, 'NoStart' => $noStart, 'NoStop' => $noStop);
#		print "\ncds $tid\t$location\n";
#		$cds->displayHeader(\*STDOUT), "\n";
		if(!defined($xcripts{$tid})) {
		    $xcripts{$tid} = [$cds, undef];
		}
		else {
		    $xcripts{$tid}->[0] = $cds;
		}
	    }
	}
    }

    return \%genes;
}

sub fuseUtrRegions {
    my($genes) = @_;
    foreach my $id (keys %$genes) {
        my $utrExons = $genes->{$id}->{'3UTRExons'};
        if(scalar(@{$utrExons}) > 1) {
            $utrExons->[0]->[2] = $utrExons->[-1]->[2];
            for(my $i = 1; $i < scalar(@{$utrExons});) {
                pop(@{$utrExons});
            }
        }
        $utrExons = $genes->{$id}->{'5UTRExons'};
        if(scalar(@{$utrExons}) > 1) {
            $utrExons->[0]->[2] = $utrExons->[-1]->[2];
            for(my $i = 1; $i < scalar(@{$utrExons});) {
                pop(@{$utrExons});
            }
        }
    }
}

sub trimStopCodon {
    my($genes) = @_;

    foreach my $id (keys %$genes) {
	my $frame =  (3 - $genes->{$id}->{'Phase'}) % 3;
	my $endPhase = (length($genes->{$id}->{'Seq'}) - $frame) % 3;
	my $lastCodon = $endPhase ? "nnn" : substr($genes->{$id}->{'Seq'},length($genes->{$id}->{'Seq'}) - 3, 3);
        if($lastCodon =~ /TGA/i || $lastCodon =~ /TAA/i || $lastCodon =~ /TAG/i) {
#       print join(" ",@{$genes{$id}->{'Exons'}->[-1]\}),"\n";
            $genes->{$id}->trimStopCodon();
        }
    }
}

sub extractIds {
    my($type, $file) = @_;
    my %ids = ();
    my $line;
    my %types = (
		 'fasta' => '^>(\S+)',
		 'gtf' => '^(\S+)\.gb\.fa',
		 'genbank' => '^LOCUS\s+([^ ]+)',
		 );
    while(defined($line = <$file>)) {
	if($line =~ /$types{$type}/) {
	    $ids{$1}++;
	}
    }
    my @res = keys %ids;
    return \@res;
}

sub revCompLocs {
    my($genes) = @_;
    my %newGenes = ();
    foreach my $gene (keys %$genes) {
	my $contig = $genes->{$gene}->{'Contig'}; 
	defined($contig) || die "Can't locate contig for $gene\n";
	$newGenes{$gene} = $genes->{$gene};
	my @exons = @{$newGenes{$gene}->{'Exons'}};
	for(my $i = 0; $i < scalar(@exons); $i++) {
	    $exons[$i]->[1] = $contig->{'Len'} - $exons[$i]->[1] + 1;
	    $exons[$i]->[2] = $contig->{'Len'} - $exons[$i]->[2] + 1;
 	}
    }
    return \%newGenes;
}

sub parseGCClasses {
    my($GCClasses)  = @_;
    my %gcClasses = ();

    my ($sep, $lastStep) = ("",0);
    my $file;
    my $class;
    foreach $sep (split(/\,/,$GCClasses)) {
        if($sep == 0) {
            next;
        }
        $gcClasses{$sep} = [$lastStep, $sep];
        $lastStep = $sep;
    }
    if(!defined($gcClasses{100})) {
        $gcClasses{100} = [$lastStep, 100];
    }
    return \%gcClasses;
}

# scramble training data for cross-validation procedure
sub scramble {
    my($type, $file, $numParts) = @_;
    my $part;
    my @IDS = @{&extractIds($type, $file)};
    my @IDS2 = ();
    my $i = 0;
    my $size = $#IDS; 
    srand;
    for ($i=0; $i<=$size; $i++)  {
	my $swap = int(rand($size+1));
	my $tmp=$IDS[$i];
	$IDS[$i] = $IDS[$swap];
	$IDS[$swap] = $tmp;
    }
    
    $i=0;
    for ($part=0; $part<$numParts; $part++)  {
	for (my $j=0; $j<($size/$numParts); $j++,$i++)
	{
	    push(@{$IDS2[$part]}, $IDS[$i]);
	}
    }
#put possible remaining few in final bin                                      
    for (;$i<=$size;$i++)  {
	push(@{$IDS2[$part - 1]}, $IDS[$i]);
    }
    return \@IDS2;
}

# Split training data for cross-validation procedure
sub  split {
    my($numSet, $IDS) = @_;
    my @trainSet = ();
    for(my $i = 0; $i < scalar(@$IDS); $i++) {
	if($i == $numSet) {
	    next;
	}
	for(my $j = 0; $j < scalar(@{$IDS->[$i]}); $j++) {
	    push(@trainSet, $IDS->[$i]->[$j]);
	}
    }
    return (\@trainSet, $IDS->[$numSet]);
}


sub markIGenics {
    my($contig, $sIgenic, $sNonIgenic, $set, $extra) = @_;
    my @array = ();
    my $i;
    my $genes = $contig->{'Features'}->{$set};
    $extra = 0 unless defined($extra);

    for($i = 0; $i <= $contig->{'Len'}; $i++) {
        push(@array, $sIgenic);
    }
    my %codRegs = ();
    my $min;
    my $max;

    while(my($id, $gene) = each (%{$genes})) {
        $min = $gene->lowestCoord() - $extra;
        $max = $gene->highestCoord() + $extra;
        
        if($min < 1 || $max < 1) {
            print STDERR "gene $id contains no exons [$max $min]\n";
            next;
        }
        
        if(defined($gene->{'gId'})) {
#            $id = $gene->{'gId'};
        }
        
        if(!defined($codRegs{$id})) {
            $codRegs{$id} = [$min, $max];
        }
        else {
            if($min < $codRegs{$id}->[0]) {
                $codRegs{$id}->[0] = $min;
            }
            if($max > $codRegs{$id}->[1]) {
                $codRegs{$id}->[1] = $max;
            }
        }
#        print STDERR "$gene->{'Id'} \"$min\" \"$max\" $contigs{$c}->{'Len'}\n";
    }
    
    while(my($gid, $range) = each (%codRegs)) {
        for($i = $range->[0]; $i <= $range->[1]; $i++) {
            $array[$i] = $sNonIgenic;
        }
    }
    return \@array;
}

sub flankingIGenicBoundaries {
    my($contig, $set) = @_;
    my $genes = $contig->{'Features'}->{$set};
    my $min = $contig->{'Len'} + 1;
    my $max = -1;

    foreach my $id (keys %$genes) {
        my $gene = $genes->{$id};
        my $exons = $gene->{'Exons'};
        
        if($gene->{'Strand'} == STRAND_FWD) {
            $min = $exons->[0]->[1];
            $max = $exons->[-1]->[2];
        }
        else {
            $min = Utils::min($min, $exons->[-1]->[2]);
            $max = Utils::max($max, $exons->[0]->[1]);
        }
    }
    return [$min, $max];
}


sub extractFlankingIGenics {
    my($contigs, $set) = @_;
    my %entries;
    foreach my $cKey (keys %$contigs) {
        my $contig = $contigs->{$cKey};
        my ($min, $max) = @{&flankingIGenicBoundaries($contig, $set)};

        if($min == 0 || $max == 0) {
            print STDERR "Warning: no genes for $contig->{'Id'}!\n";
            next;
        }
        my $fivePSeq = $contig->getSequence(0, 1, $min - 1, 0, 0, STRAND_FWD);
        my $threePSeq = $contig->getSequence(0, $max + 1, $contig->{'Len'}, 0, 0, STRAND_FWD);
        my $class = $contig->{'GCClass'};
        push(@{$entries{$class}}, [$contig->{'Id'}, '5P', $fivePSeq, $min]);
        push(@{$entries{$class}}, [$contig->{'Id'}, '3P', $threePSeq, $max]);
    }
    return \%entries;
}


sub produceAltSplices {
    my $genes = shift @_;
    my $id;
    my %altSplices = ();
    foreach $id (keys %$genes) {
        if($id =~ /(\S+)\|([^\|]+)$/) {
            push(@{$altSplices{$1}}, $2);
        }
        else {
            push(@{$altSplices{$id}}, '');
        }
    }
    return \%altSplices;
}

sub numPhaseOccurrences {
    my @occurrences = @{&phaseOccurrences(@_)};
    return $#occurrences + 1;
}

sub phaseOccurrences {
    my $seq = shift @_;
    my $phase = shift @_;
    my @substrs = @_;
    my $length = length($seq);
    my $frame =  (3 - $phase) % 3;
    my @occurrences = ();
    for(my $i = $frame; $i <= $length - 3; $i += 3) {
        my $subs = substr($seq, $i, 3);
	my $occurs = 0;
	for(my $j = 0; $j <= $#substrs; $j++) {
	    if($subs =~ /$substrs[$j]/i) {
		$occurs = 1;
		last;
	    }
	}

	if($occurs) {
#            print STDERR "$frame $i $subs\n";
	    push(@occurrences, $i);
	}
    }
    return \@occurrences;
}

sub minPhaseOccurrences {
    my $seq = shift @_;
    my @substrs = @_;
    my $minPhase;
    my $minPhaseOccurrences = length($seq);
    for(my $phase = 0; $phase < 3; $phase++) {
        my $numOccurrences = numPhaseOccurrences($seq, $phase, @substrs);
#        print "what the ", $frame,"\n";
        if($numOccurrences < $minPhaseOccurrences) {
            $minPhaseOccurrences = $numOccurrences;
            $minPhase = $phase;
        }
    }

    return ($minPhase, $minPhaseOccurrences);
}
