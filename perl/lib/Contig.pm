package Contig;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
use Gene;
use SigUtils;
use Signal;
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
	     adjustBoundsPlus
	     adjustBoundsMinus
	     getSequence
	     getSignals
	     scoreContig
	     );

use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
use enum qw(COMP=4);

sub new {
    my $type = shift;
    my %params = @_;
    my %features = ();
    
    if(defined($params{'File'})) {
        my $file = $params{'File'};
        my $line = "";

        while(defined($line = <$file>)) {
            chomp($line);

#            if($line =~ /^>?(\S+)([^\n]*)\n*([\S,\s,\n]*)$/) {
            if($line =~ /^>?(\S+)([^\n]*)\n*([\S,\n]*)$/) {
                $params{'Id'} = $1;
                $params{'Legend'} = $2;
                $params{'Seq'} = $3;
                $params{'Seq'} =~ s/\s+//g;
                last;
            }
        }
        if(!defined($line)) {
            return undef;
        }
    }

    my $self = {
	'Id' => $params{'Id'},
	'Seq' => $params{'Seq'},
      'Legend' => $params{'Legend'},
	'Beg' => defined($params{'Beg'})?$params{'Beg'}:1,
	'GCClass' => $params{'GCClass'},
	'Len' => length($params{'Seq'}),
	'Features' => \%features
    };

    $self->{'End'} = $self->{'Beg'} + $self->{'Len'} - 1;

    bless $self, $type;
    return $self;
}

sub adjustBoundsPlus {
    my($self, $beg, $end) = @_;
    my($b1,$e1);
    if ($beg < 0) {
        $b1 = 0;
    }
    else {
        $b1 = $beg;
    }
    if ($end > ($self->{'Len'} + 1))  {
        $e1 = $self->{'Len'} + 1;
    }
    else {
        $e1 = $end;
    }
    return ($b1, $e1);
}

sub adjustBoundsMinus {
    my($self, $beg, $end) = @_;
    my($b1,$e1);

    if ($end < 0) {
        $e1 = 0;
    }
    else {
        $e1 = $end;
    }
    if ($beg > ($self->{'Len'} + 1))  {
        $b1 = $self->{'Len'} + 1;
    }
    else {
        $b1 = $beg;
    }
    return ($b1,$e1);
}

sub getSequence {
    my ($self, $ups, $beg, $end, $dws, $abs, $strand) = @_;
    my ($b1, $e1) = (0,0);
    my $correction = $self->{'Beg'} - 1;
    my $seq = "";
    if(!defined($abs) || $abs == 0) {
	$correction = 0;
    }
    if(!defined($strand)) {
        $strand = STRAND_COMP;
        if($end >= $beg) {
            $strand = STRAND_FWD;
        } 
    }
#    print STDERR "Contig $self->{'Id'} : $ups, $beg, $end, $dws, $abs, $correction \n";
    my $len = $end - $beg + 1;
    if($strand == STRAND_COMP) {
	$len = $beg - $end + 1;
    }
    if(!$len) {
	return $seq;
    }	
    if($beg <= $end && $strand != STRAND_COMP)    {
        ( $b1,$e1 ) = $self->adjustBoundsPlus ($beg-$ups-$correction, $end+$dws - $correction);
        if ( not defined $b1 or not defined $e1 )       {
            die "contig $self->{'Id'}: b1 and e1 $b1,$e1\n";
        }
        elsif ( $b1 > $e1 )     {
            die "contig $self->{'Id'}: b1 and e1 $b1 > $e1";
	}
	else {
	    $seq = substr( $self->{'Seq'}, $b1-1, $e1+1-$b1 );
	}
    }
    elsif($end <= $beg && $strand != STRAND_FWD) {
	( $b1,$e1 ) = $self->adjustBoundsMinus ($beg+$ups-$correction, $end-$dws-$correction);
	if ( not defined $b1 or not defined $e1 )       {
	    die "contig $self->{'Id'}: b1 and e1 $b1, $e1";
	}
	elsif ( $e1 > $b1 )     {
	    die "contig $self->{'Id'}: e1 > b1 $e1 > $b1";
	}
	else {
	    $seq = substr( $self->{'Seq'}, $e1-1, $b1+1-$e1 );
	    $seq = &GenUtils::revComp($seq);
	}
    }
    else {
        die "Incompatible parameters in call to getSequence Contig $self->{'Id'} : $beg, $end, $strand, $self->{'Len'}\n";;
    }
    
    return $seq;
}


sub getSignals {
    my ($self, $signals, $signalPatterns) = @_;
    my $signal3bp;
    my $signal2bp;
    my $signalInd = undef;
    my $posSignal = undef;
    my $index = $self->{'Beg'};
    my $cKey = $self->{'Id'};
    if(!defined($signals->{$cKey})) {
        % { $signals->{$cKey}} = ();
    }
    while($index < $self->{'End'}) {
        if($index < $self->{'End'} - 1) {
            $signal3bp = $self->getSequence(0, $index, $index + 2, 0, 1, STRAND_FWD);
            $signalInd = undef;
            $posSignal = undef;
#    print $index - $self->{'Beg'} + 1, " ", $signal3bp, " ", $signal2bp, "\n";
            if(!($index % 100000)) {
                print STDERR ".";
            }
            if(defined($signalPatterns->{$signal3bp})) {
                $signalInd = $signalPatterns->{$signal3bp};
                $posSignal = $index;
                if($signalInd & COMP) { # complementary strand
                    $posSignal += 2;
                }
            }
            if(defined($signalInd) && defined($posSignal)) {
                if(!defined($signals->{$cKey}->{$posSignal})) {
                    @ { $signals->{$cKey}->{$posSignal}} = ();
                }
                $signals->{$cKey}->{$posSignal}->[(($signalInd & COMP) == 4)?1:0] = [$signalInd&~COMP, 0, -1];
#print "$cKey Signal ", $signalInd," ", $posSignal, " ", (($signalInd & COMP) == 4)?1:0, "\n";
            }
        }
        $signal2bp = $self->getSequence(0, $index, $index + 1, 0, 1, STRAND_FWD);
        $signalInd = undef;
        $posSignal = undef;
                                                                               
        if(defined($signalPatterns->{$signal2bp})) {
            $signalInd = $signalPatterns->{$signal2bp};
            $posSignal = $index;
            if($signalInd & COMP) { # complementary strand
                $posSignal += 1;
            }
        }

	if(defined($signalInd) && defined($posSignal)) {
            if(!defined($signals->{$cKey}->{$posSignal})) {
                @ { $signals->{$cKey}->{$posSignal}} = ();
            }
            $signals->{$cKey}->{$posSignal}->[(($signalInd & COMP) == 4)?1:0] = [$signalInd&~COMP, 0, -1];
#print "$cKey Signal ", $signalInd," ", $posSignal, " ", (($signalInd & COMP) == 4)?1:0, "\n";
        }
        $index++;
    }
}

sub scoreContig {
    my($self, $model, $gcClasses, $scoreDir) = @_;
    my $gcClass = $self->{'GCClass'};
    my $c = $self->{'Id'};
    if(!defined($gcClass)) {
	$gcClass = SigUtils::getGCClass($self->{'Seq'}, $gcClasses);
	$self->{'GCClass'} = $gcClass;
    }
    open(HANDLE, ">$scoreDir/$c.in") || die "$scoreDir/$c.in can't be opened\n";
    print HANDLE ">$c\:",STRAND_FWD,"\:", $self->{'Beg'},"\n", $self->getSequence(0, 1, length($self->{'Seq'}), 0, 1, STRAND_FWD),"\n";
    print HANDLE ">$c\:",STRAND_COMP,"\:", length($self->{'Seq'}),"\n", $self->getSequence(0, length($self->{'Seq'}), 1, 0, 1, STRAND_COMP),"\n";
    close(HANDLE);
    (-e "$model->{'CodingModel'}.$gcClass" && ! -z "$model->{'CodingModel'}.$gcClass")|| die "$model->{'CodingModel'}.$gcClass doesn't exist\n";
    my $command = "score_contigs -d $model->{'Domain'} -f $scoreDir/$c.in -M $model->{'CodingModel'}.$gcClass > $scoreDir/$c.immCoding";
    system($command);
    return $gcClass;
}
1;










