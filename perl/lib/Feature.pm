package Feature;
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Utils;
use GenUtils;
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
             dump
             getExonsWeight
             displayHeader
             displayLoc 
             displaySeq
             equalTo
             hasSameExonCoords
             Gap
             toComplement
             getFeatureSeq
             shiftExonCoords
             cropExonCoords
             genomicSeq
             );

use strict;
### Enumerates ###
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

sub new {
    my $type  = shift;
    my %params = @_;
    my $self = { };
    bless $self, $type;
    $self->initialize(%params);
    return $self;
}

sub initialize {
    my $self = shift;
    my %params = @_;
    $self->{'Id'} = $params{'Id'};
    $self->{'Seq'} = $params{'Seq'};
    my @exons = sort {&Feature::exonCompareTo($a, $b);} (@{$params{'Exons'}});
    $self->{'Exons'} = \@exons;
    $self->{'Contig'} = $params{'Contig'};
    $self->{'Strand'} = $params{'Strand'};
    $self->{'Score'} = $params{'Score'};

    if(!defined($self->{'Seq'}) || length($self->{'Seq'}) == 0) {
        $self->{'Seq'} = $self->getFeatureSeq();
    }
}

sub getExonsWeight {
    my ($self, $exons) = @_;
    my $weight = 0;
    for(my $j = 0; $j < scalar(@$exons); $j++) {
	my $exon = $exons->[$j];
	$weight += defined($exon->[3]) ? $exon->[3] : 0;
    }

    return $weight;
}


sub getFeatureSeq {
    my($self) = @_;
    if(!defined($self->{'Contig'})) {
	return $self->{'Seq'};
    }

    $self->{'Seq'} = "";
    foreach my $exon (@{$self->{'Exons'}}) {
#	print STDERR "exons ", $exon->[1], " ", $exon->[2], "\n";
        $self->{'Seq'} = $self->{'Seq'}.$self->{'Contig'}->getSequence(0,$exon->[1], $exon->[2], 0, 1, $self->{'Strand'});
#	print STDERR "length ", length( $self->{'Seq'}), "\n";
    }

    return $self->{'Seq'};
}

sub exonCompareTo {
    my($a, $b) = @_;
    my $bool = $a->[0] cmp $b->[0];
    if(!$bool) {
	return $bool;
    }
    if($a->[1] < $a->[2]) { ### positive strand
	return $a->[1] <=> $b->[1];
    }
    else {
	return -1*($a->[1] <=> $b->[1]);
    }
}

sub overlaps {
    my($self, $gene) = @_;
    my($exons_a, $exons_b) = ($self->{'Exons'},$gene->{'Exons'});
    my $bool = $self->contigId() cmp $self->contigId();
    my($min_a, $max_a, $min_b, $max_b) = (0,0,0,0);

    if($bool != 0) {
        return 0;
    }

    if($self->{'Strand'} == STRAND_FWD) {
        $min_a = $exons_a->[0]->[1];
        $max_a = $exons_a->[-1]->[2];
    }
    else {
        $min_a = $exons_a->[-1]->[2];
        $max_a = $exons_a->[0]->[1];
    }
    
    if($gene->{'Strand'} == STRAND_FWD) {
        $min_b = $exons_b->[0]->[1];
        $max_b = $exons_b->[-1]->[2];
    }
    else {
        $min_b = $exons_b->[-1]->[2];
        $max_b = $exons_b->[0]->[1];
    }
    
    return ($min_a > $min_b && $max_b > $min_a) || ($min_b > $min_a && $max_a > $min_b);
}


sub displayId {
    my($self, $file)  = @_;

    if(defined($self->{'Id'})) {
        print $file ">$self->{'Id'}\t";
        return;
    }

    (defined($self->{'tId'}) && defined($self->{'gId'})) || die "No Id available\n";

    print $file ">$self->{'gId'}|$self->{'tId'}\t";
}


sub displayHeader {
    my($self, $file, $contigs)  = @_;
    $self->displayId($file);
    $self->displayLoc($file, $self->{'Exons'}, $contigs);
}

sub getWeights4Display {
    my($self, $exons) = @_;
    my @weights = ();
    for(my $i = 0; $i < scalar(@{$exons}) ; $i++) {    
	my $exon = $exons->[$i];
	push(@weights, $exon->[3]) unless !defined($exon->[3]);
    }

    return \@weights;
}

sub displayLoc {
    my($self, $file, $exons, $contigs, $main_sep, $first_sep, $sec_sep) = @_;

    $main_sep = ";" if !defined($main_sep);
    $first_sep = "_" if !defined($first_sep);
    $sec_sep = "_" if !defined($sec_sep);
    
    my $contigBeg = (!defined($contigs) ? 1 :
		     $contigs->{$exons->[0]->[0]}->{'Beg'});
    my $i;
    my $exon;

    for($i = 0; $i < scalar(@{$exons}) ; $i++) {
	$exon = $exons->[$i];
	defined($exon->[2]) or die "$self->{'Id'} has undef exon #$i\n";
	print $file $exon->[0].$first_sep.($exon->[1] - $contigBeg + 1).$sec_sep.($exon->[2] - $contigBeg + 1);
	
	if($i < scalar(@{$exons}) - 1) {
	    print $file $main_sep;
	}
    }
    
}

sub displaySeq {
    my($self, $file) = @_;
    my $i = 0;
    if(defined($self->{'Seq'}) && length($self->{'Seq'}) > 0) {
	for($i = 0; $i < length($self->{'Seq'}); $i++) {
#      if(!($i % 80) && $i != 0) {
#          print $file "\n";
#      }
	    print $file substr($self->{'Seq'}, $i, 1);
	}
	print "\n";
    }
}

sub dump {
    my($self, $file, $contigs, $onlyHeader) = @_;
    $onlyHeader = 0 if !defined($onlyHeader);
    $self->displayHeader($file, $contigs);
    print $file "\n";
    if(!$onlyHeader) {
        $self->getFeatureSeq();
	$self->displaySeq($file);
    }
}

sub hasSameExonCoords {
    my ($self, $feature, $exonSet) = @_;
    if(scalar(@{$self->{$exonSet}}) != scalar(@{$feature->{$exonSet}})) {
        return 0;
    }
    for(my $i = 0; $i < scalar(@{$self->{$exonSet}}); $i++) {
        my $exonSelf = $self->{$exonSet}->[$i];
        my $exonFeat = $feature->{$exonSet}->[$i];
        if($exonSelf->[0] ne $exonFeat->[0] || $exonSelf->[1] != $exonFeat->[1] || $exonSelf->[2] != $exonFeat->[2]) {
            return 0;
        }
    }
    return 1;
}
    
sub equalTo {
    my ($self, $feature) = @_;
    return $self->hasSameExonCoords($feature, 'Exons');
}

### computes the gap betwee self and feature. If feature is after self, gap is positive, otherwise is negative ###
sub Gap {
    my($self, $feature) = @_;
    my($min_self, $max_self, $min_feat, $max_feat) = (0,0,0,0);
    if($self->contigId() ne $feature->contigId()) {
	return [0, 0];
    }
    if($self->{'Strand'} == STRAND_FWD) {
	if($feature->{'Strand'} == STRAND_FWD) {
	    return $feature->{'Exons'}->[0]->[1] - $self->{'Exons'}->[scalar(@{$self->{'Exons'}}) - 1]->[2];
	}
	return [$feature->{'Exons'}->[scalar(@{$feature->{'Exons'}}) - 1]->[2],  $self->{'Exons'}->[scalar(@{$self->{'Exons'}}) - 1]->[2]];
	
    }
    else {
	if($feature->{'Strand'} == STRAND_FWD) {
	    return $feature->{'Exons'}->[0]->[1] - $self->{'Exons'}->[0]->[1];
	}
	return [$feature->{'Exons'}->[scalar(@{$feature->{'Exons'}}) - 1]->[2], $self->{'Exons'}->[0]->[1]];
    }
    die "Can't happen";
}

sub toComplement {
    my($self, $lenContig) = @_;
    for(my $i = 0; $i < scalar(@{$self->{'Exons'}}); $i++) {
	$self->{'Exons'}->[$i]->[1] = $lenContig - $self->{'Exons'}->[$i]->[1];
	$self->{'Exons'}->[$i]->[2] = $lenContig - $self->{'Exons'}->[$i]->[2];
    }
    $self->{'Strand'} = ($self->{'Strand'} == STRAND_FWD) ? STRAND_COMP : STRAND_FWD;
}

sub shiftExonCoords {
    my($self, $exons, $shift, $contigLen, $newContigId) = @_;
    my @removed = (0, 0, 0);

    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
        my $exon = $exons->[$i];
#        print STDERR "shifting $exon->[1] $exon->[2] $shift\n";
        $exon->[0] = $newContigId;
        $exon->[1] = $exon->[1] + $shift;
        $exon->[2] = $exon->[2] + $shift;
       
 
        my $b = 1; my $e = 2;      
        if($self->{'Strand'} == STRAND_COMP) {
            $b = 2; $e = 1;
        }
        if($exon->[$b] < 1) {
            if($exon->[$e] < 1) {
                $removed[$b] += abs($exon->[$b] - $exon->[$e]) + 1;
                splice(@{$exons}, $i, 1);
                $i--;
                next;
            }
            else {
                $removed[$b] += 1 - $exon->[$b];
                $exon->[$b] = 1;
            }
            
        }
        if($exon->[$e] > $contigLen) {
            if($exon->[$b] > $contigLen) {
                $removed[$e] += abs($exon->[$b] - $exon->[$e]) + 1;
                splice(@{$exons}, $i, 1);
                $i--;
                next;
            }
            else {
                $removed[$e] += $exon->[$e] - $contigLen;
                $exon->[$e] = $contigLen;
            }
        }
    }

    return \@removed;
}

sub cropExonCoords {
    my($self, $exons, $lowerL, $upperL) = @_;
    my @removed = (0, 0, 0);

    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
        my $exon = $exons->[$i];
        my $b = 1; my $e = 2;      
        if($self->{'Strand'} == STRAND_COMP) {
            $b = 2; $e = 1;
        }
        
        if($exon->[$b] < $lowerL) {
            if($exon->[$e] < $lowerL) {
#		print STDERR "removing exon $exon->[$b] $lowerL\n";
                $removed[$b] += abs($exon->[$b] - $exon->[$e]) + 1;
                splice(@{$exons}, $i, 1);
                $i--;
                next;
            }
            else {
                $removed[$b] += $lowerL - $exon->[$b];
                $exon->[$b] = $lowerL;
            }
        }
        if($exon->[$e] > $upperL) {
            if($exon->[$b] > $upperL) {
#		print STDERR "removing exon $exon->[$b] $upperL\n";
                $removed[$e] += abs($exon->[$b] - $exon->[$e]) + 1;
                splice(@{$exons}, $i, 1);
                $i--;
                next;
            }
            else {
                $removed[$e] += $exon->[$e] - $upperL;
                $exon->[$e] = $upperL;
            }
        }
    }

    return \@removed;
}


sub genomicSeq {
    my($self) = @_;

    my $seq = $self->{'Contig'}->getSequence(0, $self->{'Exons'}->[0]->[1], $self->{'Exons'}->[-1]->[2], 0 , 0, $self->{'Strand'});
    return $seq;
}
1;

