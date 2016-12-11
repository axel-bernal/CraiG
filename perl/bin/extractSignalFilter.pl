#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $contigFile;
my $usage = "usage : $0 [--phase=PHASE] [--sig-weight WEIGHT] [-truncated] [-canonical] -ctg contigFile < locations > filter\n";
my $canonical = 0;
my $truncated = 0;
my $signal_weight = 1;
my $seqExon;
my $fixed_phase = -1;
&GetOptions(
    "-ctg=s" => \$contigFile,
    "-canonical" => \$canonical,
    "-phase=i" => \$fixed_phase,
    "-truncated" => \$truncated,
    "-sig-weight=f" => \$signal_weight,
    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};

for my $c (keys %contigs) {
    my @sigArray = ();
    my $i;
    my $j;
    my $gene;
    my $sigPos;
    my $donor;
    my $acceptor;
    my $clen = $contigs{$c}->{'Len'};
    for($i = 0; $i < 2; $i++) {
        @{$sigArray[$i]} = ();
        for($j = 0; $j <= $clen; $j++) {
            my @tmp = (); push(@tmp, "-");
            push(@{$sigArray[$i]}, \@tmp);
        }
    }

    my @genes = sort {&Gene::compareTo($a, $b);} (values %{ $contigs{$c}->{'Features'}->{'TRAIN'}});
    my $transcript = "";
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

	$seqExon = $gene->{'Contig'}->getSequence(0, $e_i->[1], $e_i->[1], 2, 1, $strand);	
	if(!$gene->{'NoStart'} && length($seqExon) == 3 && $seqExon =~ /atg/i) {
	    ($phase == 0) || warn "gene $gene->{'Id'} $gene->{'NoStart'} has bad start phase\n";
	    $transcript .= fillSigArray($gene, \@sigArray, "S", $e_i->[1], 0, $clen, $lastPos, $transcript, 0, $fixed_phase >= 0 ? $fixed_phase : 0, $e_i->[3], $strand);
	    $lastPos = $e_i->[1];
	}
	else {
	    $seqExon = $gene->{'Contig'}->getSequence(3, $e_i->[1], $e_i->[1], 0, 1, $strand);
	    if(length($seqExon) > 2 && !$truncated) {
		$acceptor = substr($seqExon, 1, 2);
		if($acceptor =~ /ag/i || !$canonical) {
		    $transcript .= fillSigArray($gene, \@sigArray, "A", $e_i->[1], -2, $clen, $lastPos, $transcript, 0, $fixed_phase >= 0 ? $fixed_phase : $phase, $e_i->[3], $strand);
		    $lastOffset = -2;
		    $lastPos = $e_i->[1];	    
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
			$transcript = fillSigArray($gene, \@sigArray, "D", $e_i_1->[2], 1, $clen, $lastPos, $lastOffset, $transcript, $fixed_phase >= 0 ? $fixed_phase : $phase, $e_i_1->[3], $strand);
			$lastPos = $e_i_1->[2];
		    }
		    else {
			$lastPos = -1;
			if($donor !~ /gt/i) {
			    my $compSeqExon = $gene->{'Contig'}->getSequence(3, $e_i_1->[2], $e_i_1->[2], 0, 1, ($strand == STRAND_FWD ? STRAND_COMP: STRAND_FWD));

			    warn "$gene->{'Id'} non-canonical donor $seqExon $e_i_1->[2] $compSeqExon\n";
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

                $seqExon = $gene->{'Contig'}->getSequence(3, $e_i->[1], $e_i->[1], 0, 1, $strand);
                
                if(length($seqExon) > 2) {     
                    $acceptor = substr($seqExon, 1, 2);
#                    ($acceptor =~ /ag/i) || warn "$gene->{'Id'} acceptor $seqExon $e_i->[1]\n";
                    if($acceptor =~ /ag/i || !$canonical) {
			$transcript = fillSigArray($gene, \@sigArray, "A", $e_i->[1], -2, $clen, $lastPos, 1, $transcript, $fixed_phase >= 0 ? $fixed_phase : $phase, $e_i->[3], $strand);
			$lastOffset = -2;
			$lastPos = $e_i->[1];	    
		    }
		    else {
			$lastPos = -1;
			if($acceptor !~ /ag/i) {
			    my $compSeqExon = $gene->{'Contig'}->getSequence(0, $e_i->[1], $e_i->[1], 3, 1, ($strand == STRAND_FWD ? STRAND_COMP: STRAND_FWD));
			    warn "$gene->{'Id'} non-canonical acceptor $seqExon $e_i->[1] $compSeqExon\n";
			}
		    }
		    $goodAlignment = 1;		    
		    
		}
            }
            
            $phase = ($phase + abs($e_i->[2] - $e_i->[1]) + 1) % 3;
	}

        if(!$goodAlignment && scalar(@{$exons}) > 1) {
            die "Not a good alignment for $gene->{'Id'}\n";
        }
        else {
            $e_i = $exons->[-1];
            $seqExon = $gene->{'Contig'}->getSequence(0, $e_i->[2], $e_i->[2], 3, 1, $strand);
            if(!$gene->{'NoStop'} && length($seqExon) == 4 && ($seqExon =~ /\Staa/i || $seqExon =~ /\Stga/i || $seqExon =~ /\Stag/i)) {
                ($phase == 0) || warn "gene $gene->{'Id'} $gene->{'NoStop'} has bad stop phase\n";

                $transcript .= fillSigArray($gene, \@sigArray, "T", $e_i->[2], 1, $clen, $lastPos, $lastOffset, $transcript, $fixed_phase >= 0 ? $fixed_phase : 0, $e_i->[3], $strand);
            }
            elsif(length($seqExon) > 2 && !$truncated) {
		$donor = substr($seqExon, 1, 2);
#                    ($donor =~ /gt/i || $donor =~ /gc/i) || warn "$gene->{'Id'} donor $seqExon $e_i->[2]\n";
		if($donor =~ /gt/i || !$canonical) {
		    $transcript .= fillSigArray($gene, \@sigArray, "D", $e_i->[2], 1, $clen, $lastPos, $lastOffset, $transcript, $fixed_phase >= 0 ? $fixed_phase : $phase, $e_i->[3], $strand);
		}
		elsif($donor !~ /gt/i) {
		    warn "$gene->{'Id'} non-canonical donor $seqExon $e_i->[2]\n";
		}
            }
	}
    }    
    print ">$c\t$clen\t3\n";
    GenUtils::printSigArray($sigArray[0], \*STDOUT);
    print "\/\/\n";
    GenUtils::printSigArray($sigArray[1], \*STDOUT);
}

sub fillSigArray {
    my($gene, $array, $symbol, $pos, $offsetP, $cLen, $lastPos, $lastOffset, $transcript, $phase, $weight, $strand) = @_;

    $pos = (($strand == STRAND_FWD) ? $pos + $offsetP: $cLen - ($pos - $offsetP) + 1);
    if($lastPos > 0) {
        $lastPos = (($strand == STRAND_FWD) ? $lastPos + $lastOffset: $cLen - ($lastPos - $lastOffset) + 1);
    }

    my $key = "$lastPos $phase";
    my $this_weight = defined($weight) ? $weight : $signal_weight;
    GenUtils::fillSigFastaArray($gene, $this_weight, $array, $symbol, $pos, $key, $strand);
}
