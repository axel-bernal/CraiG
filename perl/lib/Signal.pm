package Signal;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
use SigUtils;
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
	     initSeq
	     getAssocExon
	     );

use enum qw(START STOP DONOR ACEPTOR);
use enum qw(COMP=4);     
use enum qw(STRAND_FWD STRAND_COMP);
use strict;

sub new {
    my $type = shift;
    my %params = @_;
    my $self = {
	'GeneId' => $params{'Geneid'},
	'Pos' => $params{'Pos'}, # relative to the beginning of the contig. format is "contig"._."pos"
	'Seq' => $params{'Seq'},
	'Type' => $params{'Type'},
	'ExonLen' => $params{'ExonLen'},
	'ExonGCC' => $params{'ExonGCC'},
	'Class' => $params{'Class'},
	'PredClass' => $params{'PredClass'},
	'Score' => $params{'Score'},
	'Strand' => $params{'Strand'},
	'ForTraining' => $params{'ForTraining'},
	'State' => $params{'State'},
    };
    bless $self, $type;
    return $self;
}

sub initSeq {
    my($self, $contigs, $window) = @_;
    $window =~ /(\d+)\:(\d+)/ || die "Bad format $window\n";
    my $beforeSignal = $1;
    my $afterSignal = $2;
    $self->{'Pos'} =~ /(\S+)\_(\d+)/;
    my $pSCoord = $2;
    my $pSContig = $1;
    
    my $signalSeq = undef;
    if($self->{'Strand'} == STRAND_COMP) { # complementary strand
	if($pSCoord + $beforeSignal > $contigs->{$pSContig}->{'End'} || $pSCoord - $afterSignal < 1) {
	    $self->{'ForTraining'} = 0;
	    return;
	}
	$signalSeq = $contigs->{$pSContig}->getSequence($beforeSignal, $pSCoord, $pSCoord - 1, $afterSignal - 2, 1, $self->{'Strand'});
    }
    else {
	if($pSCoord - $beforeSignal < 1 || $pSCoord + $afterSignal > $contigs->{$pSContig}->{'End'}) {
	    $self->{'ForTraining'} = 0;
	    return;
	}
	$signalSeq = $contigs->{$pSContig}->getSequence($beforeSignal, $pSCoord, $pSCoord + 1, $afterSignal - 2, 1, $self->{'Strand'});
    }
    # computing additional signal information
    $self->{'Seq'} = $signalSeq;
#    print join(" ", %$self),"\n";
}

sub getAssocExon {
    my($self, $numExons, $maxLen) = @_;
    $self->{'Pos'} =~ /(\S+)\_(\d+)/;
    my $pSCoord = $2;
    my $pSContig = $1;
    my $exon = [$pSContig, $pSCoord, $pSCoord];
    my $dir = 1;
    # changing shift and step according to the strand and signal info
    if($self->{'Strand'} == STRAND_COMP) {
#	print "hola\n";
	if(($self->{'Type'} == START) || ($self->{'Type'} == ACEPTOR)) {
	    $dir = -1;
	}
    }
    else { # ends in a stop, shift will always be zero
	if(($self->{'Type'} == DONOR) || ($self->{'Type'} == STOP)) {
            $dir = -1;
	}
    }
#    print "Dir $dir Strand $self->{'Type'} MaxLen $maxLen ", $self->{'Type'} & ~COMP,"\n";
    my $len = $dir*int(rand($maxLen));
    $self->{'Len'} = abs($len + $dir*1);
    $exon->[2] += $len;
    return $exon;
}
1;
