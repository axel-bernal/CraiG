#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $contigFile = "";
my $locations = "";
my $canonical = 0;
my $normCounts = 0;

my $usage = "usage : $0 [-canonical] [-norm_counts] -ctg contigFile -locs weighted_locations -canonical.\n\tWrites all gernerated filters in filter files *.stscore and *.state\n";
my $seqExon;

&GetOptions(
    "-locs=s" => \$locations,
    "-ctg=s" => \$contigFile,
    "-canonical" => \$canonical,
    "-norm_counts" => \$normCounts,
    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;
(defined($locations) && length($locations) > 0) || die "Locations File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my $one_gene = 1;
my $no_alt = 1;
open(LOCS, "<$locations.locs") || die "file $locations not found\n";
my $genes = &GenUtils::readWeightedFastaFormat(\*LOCS, 0, 'TRAIN', \%contigs);
close(LOCS);

if($normCounts) {
    open(NORM_SCORE, ">$locations.norm.stscores") || die "Cannot open $locations.norm.stscores file\n";
}

open(STATE, ">$locations.state") || die "Cannot open $locations.state file\n";
open(SCORE, ">$locations.stscores") || die "Cannot open $locations.stscores file\n";
open(SIGNAL,">$locations.addsignal") || die "Cannot open $locations.addsignal file\n";

my $maxCounts = 0;

for my $c (keys %contigs) {

    my @stArray = ();
    my @sigArray = ();
    my @signals = ();
    my @countArray = ();
    my @maxArray = ();
    my $i;
    my $j;
    my $gene;
    my $sigPos;
    my $donor;
    my $acceptor;
    $maxCounts = 0;
    
    for($i = 0; $i < 2; $i++) {
        @{$sigArray[$i]} = ();
        @{$stArray[$i]} = ();
        @{$signals[$i]} = ();
	@{$countArray[$i]} = ();
	@{$maxArray[$i]} = ();
	
        for($j = 0; $j <= $contigs{$c}->{'Len'}; $j++) {
            my @tmp = (); push(@tmp, "-");
            push(@{$sigArray[$i]}, \@tmp);
            push(@{$stArray[$i]}, ".");
	    push(@{$countArray[$i]}, 0);
	    push(@{$maxArray[$i]}, 0);
        }
    }
    
    my $cgenes = $contigs{$c}->{'Features'}->{'TRAIN'};
    for my $mcid (keys %$cgenes) {
	$cgenes->{$mcid}->fixAnnot4Alignment();
    }

    my @genes = sort {&Gene::compAlignQualityTo($a, $b);} (values %{ $cgenes});

    for $gene (@genes) {

        my $goodAlignment = 0;
        my $strand = $gene->{'Strand'};
        my $phase = $gene->{'Phase'};
	my $exons = $gene->{'Exons'};
	if(scalar(@{$exons}) == 0) {
	    next;
	}

        for($i = 0; $i < scalar(@{$exons}); $i++) {
	    my $fillChar = '.';
	    my $upsChar = '';
	    my $dsChar = '';
	    my $e_i = $exons->[$i];
	    print STDERR "receiving $gene->{'Id'} $e_i->[1] $e_i->[2]\n";
	    my $phase = $e_i->[4];
	    my $numStops = $e_i->[5];

	    if($i == 0 && !$gene->{'NoStart'} && $e_i->[6]) {
		($phase == 0) || warn "gene $gene->{'Id'} $gene->{'NoStart'} has bad start phase\n";
		$upsChar = '.';
		$fillChar = 'S';
		if(scalar(@$exons) > 1) {
		    $fillChar = 'I';
		}
	    }
	    elsif($e_i->[6] || !$canonical) {
		$fillChar = 'E';
		if($e_i->[6]) {
		    $upsChar = 'i';
		}
	    }
	    else {
		print STDERR "$e_i->[1] not an start/acceptor\n";
		next;
	    }
	    
	    my $beg;
	    my $end;

	    if($i > 0) { #fill intron
		my $e_i_1 = $exons->[$i - 1];
		$beg = (($strand == STRAND_FWD) ?
			$e_i_1->[2] + 1:
			$contigs{$c}->{'Len'} - $e_i_1->[2] + 2);
		$end = (($strand == STRAND_FWD) ?
			$e_i->[1] - 1: 
			$contigs{$c}->{'Len'} - $e_i->[1]);

		fillStArray($gene, 0, $countArray[$strand], $maxArray[$strand], $stArray[$strand], 'i', $beg, $end, $contigs{$c}->{'Len'}, 0);
	    }
	    
            if($i == scalar(@$exons) - 1 && !$gene->{'NoStop'} && $e_i->[7]) {
		((($phase + abs($e_i->[2] - $e_i->[1]) + 1) % 3) == 0)
		    || warn "gene $gene->{'Id'} $gene->{'NoStop'} has bad stop phase\n";
		$fillChar = 'S';
		$dsChar = '.';
		if(scalar(@$exons) > 1) {
		    $fillChar = 'L';
		}
	    }
	    elsif($e_i->[7] || !$canonical) {
#                    ($donor =~ /gt/i || $donor =~ /gc/i) || warn "$gene->{'Id'} donor $e_i->[2]\n";
		if($e_i->[7]) {
		    $dsChar = 'i';
		}
	    }
	    else {
		print STDERR "$e_i->[2] not a donor\n";
		next;
	    }

#	    if($upsChar eq '.' && $dsChar eq '.' && $numStops) {
#		print STDERR "$e_i->[1] $e_i->[2] is malformed\n";
#		next;
#	    }

	    $beg = (($strand == STRAND_FWD) ?
		    $e_i->[1]:
		    $contigs{$c}->{'Len'} - $e_i->[1] + 1);
	    $end = (($strand == STRAND_FWD) ? 
		    $e_i->[2]:
		    $contigs{$c}->{'Len'} - $e_i->[2] + 1);
	    
	    if($beg > 1 && $upsChar ne '') {
		push(@{$signals[$strand]}, [$beg - 1, $upsChar]);
	    }
	 
	    if($end < $contigs{$c}->{'Len'} && $dsChar ne '') {
		push(@{$signals[$strand]}, [$end + 1, $dsChar]);
	    }   

#	    fillStArray($gene, $e_i->[3], $countArray[$strand], $maxArray[$strand], $stArray[$strand], $fillChar, $beg, $end, $contigs{$c}->{'Len'}, $phase);
 	    fillStArray($gene, $e_i->[3], $countArray[$strand], $maxArray[$strand], $stArray[$strand], 'E', $beg, $end, $contigs{$c}->{'Len'}, $phase);           
	    $phase = ($phase + abs($e_i->[2] - $e_i->[1]) + 1) % 3;
	    
	}
       
    }
    
    for($i = 0; $i < 2; $i++) {
        for($j = 0; $j < scalar(@{$signals[$i]}); $j++) {
	    my $sig = $signals[$i]->[$j];
	    $stArray[$i]->[$sig->[0]] = $sig->[1];
	}
    }
    
# eliminate segments with multiple letters, when possible, this will be an option 
# which wont be used if the user wants information regarding alternative splicing
    print SIGNAL ">$c\t$contigs{$c}->{'Len'}\t3\n";
    for($i = 0; $i < 2; $i++) {
    	for($j = 0; $j <= $contigs{$c}->{'Len'}; $j++) {
    	    if($stArray[$i]->[$j] =~ /[I,S,L,E,J,T,M,A,-,_]/) {
    		my $c1 = $j;
    		my $c2 = $j + 2;
    		if($i == 1) {
    		    $c1 = $contigs{$c}->{'Len'} - $c1 + 1;
    		    $c2 = $contigs{$c}->{'Len'} - $c2 + 1;
    		}
    		my $stop = $contigs{$c}->getSequence(0, $c1, $c2, 0, 1, $i);
    		if($stop =~ /tag/i || $stop =~ /tga/i || $stop =~ /taa/i) {
		    GenUtils::fillSigFastaArray($gene, \@sigArray, 'T', $j, "-1 0", $i);
    		}
    	    }
    	}

	GenUtils::printSigArray($sigArray[$i], \*SIGNAL);
	if($i == 0) {
	    print SIGNAL "\/\/\n";
	}
    }
    
    print STATE ">$c\n";
    GenUtils::printStateArray($stArray[0], \*STATE);
    print STATE "\/\/\n";
    GenUtils::printStateArray($stArray[1], \*STATE);

    print SCORE ">$c\t$contigs{$c}->{'Len'}\t2\n";
    printAllFilters($countArray[0], $maxArray[0], \*SCORE);
    print SCORE "\/\/\n";
    printAllFilters($countArray[1], $maxArray[1], \*SCORE);

    # normalize counts
    if($normCounts) {
	for($i = 0; $i < 2; $i++) {
	    for($j = 0; $j <= $contigs{$c}->{'Len'}; $j++) {
		if($maxCounts) {
		    $countArray[$i]->[$j] = $countArray[$i]->[$j]*100/$maxCounts;
		}
		else {
		    $countArray[$i]->[$j] = 0;
		}
	    }
	}
	print NORM_SCORE ">$c\t$contigs{$c}->{'Len'}\t2\n";
	printAllFilters($countArray[0], $maxArray[0], \*NORM_SCORE);
	print NORM_SCORE "\/\/\n";
	printAllFilters($countArray[1], $maxArray[1], \*NORM_SCORE);
	
    }

}

close(STATE);
close(SCORE);
close(SIGNAL);

if($normCounts) {
    close(NORM_SCORE);
}

sub printAllFilters {
    my ($carray_ref, $marray_ref, $file) = @_;

    (scalar(@{$carray_ref}) > 0 && scalar(@{$marray_ref}) > 0) || die "invalid array ref\n";
    
    my $pos;
    for(my $i = 1; $i < scalar(@{$carray_ref}); $i++)  {	
        if($carray_ref->[$i] || $marray_ref->[$i]) {
            print $file "$i\t", (sprintf "%.2f", $carray_ref->[$i]),"\t", (sprintf "%.5f", $marray_ref->[$i]), "\n";
	}
    }
}

sub fillStArray {
    my($gene, $weight, $countArray, $maxArray, $array, $fillChar, $b, $e, $cLen, $phase) = @_;
    my %intron_map = ('I', 'J', 'S', 'T', 'L', 'M', 'E','A');    
    print STDERR "filling $gene->{'Id'} $b $e $phase $fillChar\n";

    my @begins = ();
    push(@begins, [$b, $phase]);
    my $pos;
    if($fillChar eq 'i') {
	$pos = $b;
	for( ; $pos <= $e; $pos++) {
	    if($array->[$pos] eq '.') {
		$array->[$pos] = $fillChar;
	    }
	}	
	return;
    }
    elsif($fillChar =~ /[I,S,L,E]/) {
	my $qpos = $b;
	for( ; $qpos <= $e; $qpos++) {
	    if($array->[$qpos] =~ /[I,S,L,E,J,T,M,A]/) {
		last;
	    }
	}

	if($qpos <= $e) {
	    $begins[-1]->[1] = (($b - $qpos) % 3 + 3) % 3;
	}

	$pos = $b;
	my @segment = ();	
	for( ; $pos <= $e; $pos++) {
	    my $opChar = '-';
	    if($array->[$pos] =~ /i/) {
		$opChar = '_';
	    }

	    push(@segment, $opChar);

	}

	$phase = $begins[-1]->[1];
	$pos = $begins[-1]->[0] + (3 - $phase) % 3; 

	for( ; $pos <= $e; $pos += 3) {
	    my $opChar = $fillChar;
	    if($array->[$pos] =~ /i/) {
		$opChar = $intron_map{$fillChar};
	    }
	    elsif($array->[$pos] =~ /[\-,\_]/) {
		my $qpos = $pos - 2;
		for( ; $qpos <= $e; $qpos++) {
		    if($qpos >= $b && $array->[$qpos] !~ /[i,\_]/) {
			last;
		    }
		}

		$pos = $qpos;
		push(@begins, [$pos, $phase]);
		$qpos = $pos;
		for( ; $qpos <= $e; $qpos++) {
		    if($array->[$qpos] =~ /[I,S,L,E,J,T,M,A]/) {
			last;
		    }
		}
		
		if($qpos <= $e) {
		    $begins[-1]->[1] = (($pos - $qpos) % 3 + 3) % 3;
		}

		$phase = $begins[-1]->[1];
		$pos = $begins[-1]->[0] + (3 - $phase) % 3; 
	    }
	    $segment[$pos - $b] = $opChar;
	}

	$pos = $b;
	for( ; $pos <= $e; $pos++) {
	    my $nletter = $segment[$pos - $b];
	    if($array->[$pos] =~ /\./) {
		$array->[$pos] = $nletter;
	    }
	    elsif(scalar(@begins) == 1 && $array->[$pos] =~ /i/) {
		$array->[$pos] = $nletter;
	    }
	    elsif($array->[$pos] =~ /[I,S,L,E,J,T,M,A]/) {
		;
	    }
	}
    }

#    print STDERR ">\n";
#    GenUtils::printStateArray($array, \*STDERR);
#    print STDERR "\n";

    my $j = 0;
    $pos = $begins[$j]->[0] + (3 - $begins[$j]->[1]) % 3;     

    for( ; $pos <= $e; $pos += 3) {
	if($j < $#begins && $pos >= $begins[$j+1]->[0]) {
	    $j++;
	    $pos = $begins[$j]->[0] + (3 - $begins[$j]->[1]) % 3; 
	}

	if($maxArray->[$pos] < $weight) {
	    $maxArray->[$pos] = $weight;
	}
	
	$countArray->[$pos]++;
	if($countArray->[$pos] > $maxCounts) {
	    $maxCounts = $countArray->[$pos];
	}	
    }        
}
