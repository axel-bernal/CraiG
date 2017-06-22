package Organism;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
require Exporter;
use Utils;
use Gene;

@ISA = ( Exporter );
@EXPORT= qw (
             loadContigs
             loadGenes
             tieContigs
             getGCClasses
             getContigSignals
             );

use strict;

### Enumerates ###
use enum qw(START STOP DONOR ACEPTOR);
use enum qw(STRAND_FWD STRAND_COMP);
use enum qw(COMP=4);

sub new  {
    my $type = shift;
    my %params = @_;
    my %contigs = ();
    my %features = ();
    my %promoters = ();
    my $self = {
	'Contigs' => \%contigs,
	'Features' => \%features,
	'Name' => $params{'Name'},
	'Domain' => $params{'Domain'},
	# last four signals are in the complementary strand
	'SignalPatterns' => {'ATG' => START, # start cdn -> index in 'Signals' array
			     'TGA' => STOP, # stop cdn
			     'TAA' => STOP,
			     'TAG' => STOP,
			     'GT' => DONOR, # donor signal
			     'AG' => ACEPTOR, # aceptor signal
			     'CAT' => START|COMP, # start cdn -
			     'TCA' => STOP|COMP, # stop cdn -
			     'TTA' => STOP|COMP,
			     'CTA' => STOP|COMP,
			     'AC' => DONOR|COMP, # donor cdn -
			     'CT' => ACEPTOR|COMP, # aceptor cdn -
			 },	
    };
    bless $self, $type;
}


sub loadContigs {
    my ($self, $contigFile) = @_;
    my $seq;
    my $id;
    my $bad;
    my $line = <$contigFile>;
    while (defined($line)) {
	my $beg = 1;
#	if($line =~ /^>(\S+)\:(\d+)/) {
#	    $id = $1; $beg = $2;
#	}
#	els
	if($line =~ /^>(\S+)/) {
	    $id = $1;
	}
	else {
	    $line = <$contigFile>; 
	    next;
	}
      chomp($id);
      $seq = "";
      while (defined($line = <$contigFile>) && ($line !~ /^>/)) {
          $seq .= $line;
      }
      $seq =~ s/\s//g;
      $seq =~ s/U/T/g;
      $seq =~ s/u/t/g;
#      if ($seq =~ /([^ACGTMRWSYKBDHVNacgtmrwsykbdhvn])/)
#      {
#          $bad = "";
#          while ($seq =~ s/([^ACGTMRWSYKBDHVNacgtmrwsykbdhvn])//)
#          {
#              $bad .= "$1 ";
#          }
#          print STDERR "not processing $id due to bad characters: $bad\n";
#      }
      $self->{'Contigs'}->{$id} = Contig->new('Id' => $id, 'Seq' => $seq, 'Beg' => $beg);
    }
    return $self->{'Contigs'};
}

sub loadGenes {
    my ($self, $file, $type, $includeSeq, $fmt) = @_;

    $type = 'TRAIN' if !defined($type);
    $includeSeq = 0 unless defined($includeSeq);
    if(!defined($fmt)) { $fmt = 'fasta';}
    if($fmt eq "gtf") {
        return($self->{'Features'}->{$type} = &GenUtils::readGFFFormat($file, $type, $self->{'Contigs'}));
    }
    elsif($fmt eq "fasta")  {
        return($self->{'Features'}->{$type} =  &GenUtils::readLocsFormat($file, $includeSeq, $type, $self->{'Contigs'}));
    }
    elsif($fmt eq "genbank") {
        return($self->{'Features'}->{$type} =  &GenUtils::readGenBankFormat(0, $file, $type, $includeSeq, $self->{'Contigs'}));
    }
    die "$fmt unknown\n";
}

sub getGCClasses {
    my ($self, $type, $gcClasses) = @_;
    my $len;
    my $location;
    my $id;
    my $contig;
    my $idClass;
    my %classes = ();
    defined($self->{'Features'}->{$type}) || return \%classes;
    my $cKey;
    
    foreach $cKey (keys %{$self->{'Contigs'}}) {
	$contig = $self->{'Contigs'}->{$cKey};
#	$idClass = SigUtils::getSeqClassBySignal('G|C', 1, $contig->{'Seq'}, $gcClasses);
	$idClass = SigUtils::getGCClass($contig->{'Seq'}, $gcClasses);
	$contig->{'GCClass'} = $idClass;
#	print "$contig->{'Id'} $contig->{'GCClass'}\n";
    }
    
    my %locs = % { $self->{'Features'}->{$type}};
    foreach $id (keys %locs) {
        $contig = $locs{$id}->{'Contig'};
        defined $contig or die $id;
#        print "$id $contig->{'Id'} $contig->{'GCClass'}\n";
        push(@ { $classes{$contig->{'GCClass'}}}, $id);
    }
    return \%classes;
}


sub getTranscriptGCClasses {
    my ($self, $type, $gcClasses) = @_;
    my $id;
    my $idClass;
    my %classes = ();
    defined($self->{'Features'}->{$type}) || return \%classes;
    my %gcClasses = ();
    my %locs = % { $self->{'Features'}->{$type}};

    foreach $id (keys %locs) {
	my $transcript = $locs{$id};
        my $c = $transcript->{'Contig'};
        defined $c or die $id;
	my $seq = $c->getSequence(0, $transcript->{'Exons'}->[0]->[1],
				  $transcript->{'Exons'}->[-1]->[2], 0, 0,  $transcript->{'Strand'});

	$idClass = SigUtils::getGCClass($seq, $gcClasses);

#        print "$id $contig->{'Id'} $contig->{'GCClass'}\n";
        $gcClasses{$id} = $idClass;
    }

    return \%gcClasses;
}

sub getContigSignals {
    my($self) = @_;
    my %signals = ();
    my @contigKeys = sort {$a cmp $b} keys (%{$self->{'Contigs'}});
    my $i;
    my $cKey;

    foreach $cKey (@contigKeys) {
        my $contig = $self->{'Contigs'}->{$cKey};
        $contig->getSignals(\%signals, $self->{'SignalPatterns'});
    }
    return \%signals;
}

1;
