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
my $truncated = 0;
my $usage = "usage : $0 [-truncated] [-canonical] [-norm_counts] -ctg contigFile -locs weighted_locations\n-canonical only considers canonical splice sites.\n\tWrites all generated filters in filter files locs.sigscore and locs.signal\n";
my $seqExon;

&GetOptions(
    "-locs=s" => \$locations,
    "-ctg=s" => \$contigFile,
    "-canonical" => \$canonical,
    "-norm_counts" => \$normCounts,
    "-truncated" => \$truncated,
    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;
(defined($locations) && length($locations) > 0) || die "Locations File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

open(LOCS, "<$locations.locs") || die "file $locations.locs not found\n";
my $genes = &GenUtils::readLocsFormat(\*LOCS, 0, 'TRAIN', \%contigs);
close(LOCS);

if($normCounts) {
    open(NORM_SCORE, ">$locations.norm.sigscores") || die "Cannot open $locations.norm.sigscores file\n";
}
open(SIGNAL, ">$locations.signal") || die "Cannot open $locations.signal file\n";
open(SCORE, ">$locations.sigscores") || die "Cannot open $locations.sigscores file\n";

my $maxCounts = 0;

for my $c (keys %contigs) {

    my @sigArray = ();
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
	@{$countArray[$i]} = ();
	@{$maxArray[$i]} = ();
	
        for($j = 0; $j <= $contigs{$c}->{'Len'}; $j++) {
            my @tmp = (); push(@tmp, "-");
            push(@{$sigArray[$i]}, \@tmp);
	    push(@{$countArray[$i]}, 0);
	    push(@{$maxArray[$i]}, 0);
        }
    }

    my @genes = sort {&Gene::compareTo($a, $b);} (values %{ $contigs{$c}->{'Features'}->{'TRAIN'}});

    for $gene (@genes) {

        my $goodAlignment = 0;
        my $strand = $gene->{'Strand'};
        my $lastOffset = 0;
	my $lastPos = -1;
        my $phase = $gene->{'Phase'};
	my $exons = $gene->{'Exons'};

	if(scalar(@{$exons}) == 0) {
	    next;
	}
	
	my $e_i = $exons->[0];
	my $this_seq = $gene->{'Contig'}->getSequence(0, $e_i->[1], $e_i->[2], 0, 1, $strand);
	
	my $numStops = GenUtils::numPhaseOccurrences($this_seq, $phase, "taa", "tga", "tag");
	if($numStops > 0) {
	    print STDERR "$gene->{'Contig'}->{'Id'} $gene->{'Id'} $e_i->[1] $e_i->[2] with $numStops stops in phase $phase\t";
	    ($phase, $numStops) = GenUtils::minPhaseOccurrences($this_seq, "taa", "tga", "tag");

	    if(!$numStops) {
		$gene->{'Phase'} = $phase;
		print STDERR "\tfixed with phase $phase\n";
	    }
	    else { print STDERR "\tnew phase $phase gives $numStops stops\n";}
	}

	if($truncated) {
	    $phase = 3;
	}

	$seqExon = $gene->{'Contig'}->getSequence(0, $e_i->[1], $e_i->[1], 2, 1, $strand);	
	if(!$gene->{'NoStart'} && length($seqExon) == 3 && $seqExon =~ /atg/i) {
	    ($phase == 0) || warn "gene $gene->{'Id'} $gene->{'NoStart'} has bad start phase\n";
	    if(GenUtils::fillSigArray($gene, $e_i->[3], \@countArray, \@maxArray, \@sigArray, "S", $e_i->[1], 0, $contigs{$c}->{'Len'}, $lastPos, 0, 0, \$maxCounts, $strand)) {
		$lastPos = $e_i->[1];
	    }
	}
	else {
	    $seqExon = $gene->{'Contig'}->getSequence(3, $e_i->[1], $e_i->[1], 0, 1, $strand);
	    if(length($seqExon) > 2 && !$truncated) {
		$acceptor = substr($seqExon, 1, 2);
		if($acceptor =~ /ag/i || !$canonical) {
		    if(GenUtils::fillSigArray($gene, $e_i->[3], \@countArray, \@maxArray, \@sigArray, "A", $e_i->[1], -2, $contigs{$c}->{'Len'}, $lastPos, 0, $phase, \$maxCounts, $strand)) {
			$lastOffset = -2;
			$lastPos = $e_i->[1];	    
		    }
		}
		elsif($acceptor !~ /ag/i) {
		    warn "$gene->{'Id'} non-canonical acceptor $seqExon $e_i->[1]\n";
		}
	    }
	}
	
        for(my $i = 0; $i < scalar(@{$exons}); $i++) {
            $e_i = $exons->[$i];

            if($i > 0) {
                my $e_i_1 = $exons->[$i - 1];
                $seqExon = $gene->{'Contig'}->getSequence(0, $e_i_1->[2], $e_i_1->[2], 3, 1, $strand);

                if(length($seqExon) > 2) {
                    $donor = substr($seqExon, 1, 2);
#                    ($donor =~ /gt/i || $donor =~ /gc/i) || warn "$gene->{'Id'} donor $seqExon $e_i_1->[2]\n";
                    if($donor =~ /gt/i || !$canonical) {
			if(GenUtils::fillSigArray($gene, $e_i_1->[3], \@countArray, \@maxArray, \@sigArray, "D", $e_i_1->[2], 1, $contigs{$c}->{'Len'}, $lastPos, $lastOffset, $phase, \$maxCounts, $strand)) {
			    $lastPos = $e_i_1->[2];
			}
		    }
		    else {
			$lastPos = -1;
			if($donor !~ /gt/i) {
			    warn "$gene->{'Id'} non-canonical donor $seqExon $e_i_1->[2]\n";
			}
		    }
		    $goodAlignment = 1;
		}

		$this_seq = $gene->{'Contig'}->getSequence(0, $e_i->[1], $e_i->[2], 0, 1, $strand);
		
		my $numStops = GenUtils::numPhaseOccurrences($this_seq, $phase, "taa", "tga", "tag");
		if($numStops > 0) {
		    print STDERR "$gene->{'Contig'}->{'Id'} $gene->{'Id'} $e_i->[1] $e_i->[2] with $numStops stops in phase $phase\t";
		($phase, $numStops) = GenUtils::minPhaseOccurrences($this_seq, "taa", "tga", "tag");

		    if(!$numStops) {
			print STDERR "\tfixed with phase $phase\n";
		    }
		    else { print STDERR "\tnew phase $phase gives $numStops stops\n";}
		}

		if($truncated) {
		    $phase = 3;
		}

                $seqExon = $gene->{'Contig'}->getSequence(3, $e_i->[1], $e_i->[1], 0, 1, $strand);
                
                if(length($seqExon) > 2) {     
                    $acceptor = substr($seqExon, 1, 2);
#                    ($acceptor =~ /ag/i) || warn "$gene->{'Id'} acceptor $seqExon $e_i->[1]\n";
                    if($acceptor =~ /ag/i || !$canonical) {
			if(GenUtils::fillSigArray($gene, $e_i->[3], \@countArray, \@maxArray, \@sigArray, "A", $e_i->[1], -2, $contigs{$c}->{'Len'}, $lastPos, 1, $phase, \$maxCounts,  $strand)) {
			    $lastOffset = -2;
			    $lastPos = $e_i->[1];	    
			}
		    }
		    else {
			$lastPos = -1;
			if($acceptor !~ /ag/i) {
			    warn "$gene->{'Id'} non-canonical acceptor $seqExon $e_i->[1]\n";
			}
		    }
		    $goodAlignment = 1;		    
		    
		}
            }
            
            $phase = ($phase + abs($e_i->[2] - $e_i->[1]) + 1) % 3;

	    if($truncated) {
		$phase = 3;
	    }
        }

        if(!$goodAlignment && scalar(@{$exons}) > 1) {
            die "Not a good alignment for $gene->{'Id'}\n";
        }
        else {
            $e_i = $exons->[-1];
            $seqExon = $gene->{'Contig'}->getSequence(0, $e_i->[2], $e_i->[2], 3, 1, $strand);
            if(!$gene->{'NoStop'} && length($seqExon) == 4 && ($seqExon =~ /\Staa/i || $seqExon =~ /\Stga/i || $seqExon =~ /\Stag/i)) {
                ($phase == 0) || warn "gene $gene->{'Id'} $gene->{'NoStop'} has bad stop phase\n";

                GenUtils::fillSigArray($gene, $e_i->[3], \@countArray, \@maxArray, \@sigArray, "T", $e_i->[2], 1, $contigs{$c}->{'Len'}, $lastPos, $lastOffset, 0, \$maxCounts, $strand);
            }
            elsif(length($seqExon) > 2 && !$truncated) {
		$donor = substr($seqExon, 1, 2);
#                    ($donor =~ /gt/i || $donor =~ /gc/i) || warn "$gene->{'Id'} donor $seqExon $e_i->[2]\n";
		if($donor =~ /gt/i || !$canonical) {
		    GenUtils::fillSigArray($gene, $e_i->[3], \@countArray, \@maxArray, \@sigArray, "D", $e_i->[2], 1, $contigs{$c}->{'Len'}, $lastPos, $lastOffset, $phase, \$maxCounts, $strand);
		}
		elsif($donor !~ /gt/i) {
		    warn "$gene->{'Id'} non-canonical donor $seqExon $e_i->[2]\n";
		}
            }
        }
    }
    
    print SIGNAL ">$c\t$contigs{$c}->{'Len'}\t3\n";
    GenUtils::printSigArray($sigArray[0], \*SIGNAL);
    print SIGNAL "\/\/\n";
    GenUtils::printSigArray($sigArray[1], \*SIGNAL);

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

close(SIGNAL);
close(SCORE);

if($normCounts) {
    close(NORM_SCORE);
}

sub printAllFilters {
    my ($carray_ref, $marray_ref, $file) = @_;

    (scalar(@{$carray_ref}) > 0 && scalar(@{$marray_ref}) > 0) || die "invalid array ref\n";
    
    my $pos;
    for(my $i = 1; $i < scalar(@{$carray_ref}); $i++)  {	
        if($carray_ref->[$i] || $marray_ref->[$i]) {
            print $file "$i\t", (sprintf "%.2f", $carray_ref->[$i]),"\t", (sprintf "%.3f", $marray_ref->[$i]),"\n";
	}
    }
}

