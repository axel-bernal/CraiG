package GeneModelProps;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
use GenUtils;
use Utils;
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
	     syncInfo
	     storeProps
	     retrieveProps
	     );

use strict;

sub new {
    my ($type) = shift;
    my %params = @_;
    if(scalar(keys %params) == 0) {
	my $self = {
	    'Name' => "",
	};
 	bless $self, $type;
	return $self;
    }

    my @signalWindows = @{$params{'SignalWindows'}};
    my @signalModels = @{$params{'SignalModels'}};
    my @signalFiles = @{$params{'SignalFiles'}};
    my @geoMeans = @{$params{'IGenicGeoMeans'}};
    my $self = {
	'Name' => $params{'Name'},
	'Path' => $params{'Path'},
	'Org' => $params{'Org'},
	'Domain' => $params{'Domain'},
	'PromoterWindow' => $params{'PromoterWindow'},
	'SignalWindows' => \@signalWindows,
	'SignalModels' => \@signalModels,
	'SignalFiles' => \@signalFiles,
        'IGenicGeoMeans' => \@geoMeans,
	'HistExonNum' => $params{'HistExonNum'},
	'HistExonLen' => $params{'HistExonLen'},
	'HistIntronLen' => $params{'HistIntronLen'},
	'HistIntergenicLen' => $params{'HistIntergenicLen'},
	'CodingModel' => $params{'CodingModel'},
	'IntronModel' => $params{'IntronModel'},
	'IntronModel' => $params{'IntronModel'},
	'IntergenicModel' => $params{'IntergenicModel'},
	'GCProp' => $params{'GCProp'},
	'GHMMParams' => $params{'GHMMParams'},
	'SigScoreBinLen' => $params{'SigScoreBinLen'},
	'CodScoreBinLen' => $params{'CodScoreBinLen'},
	'IStScoreBinLen' => $params{'IStScoreBinLen'},
	'HistIntergenicBinLen' => $params{'HistIntergenicBinLen'},
	'HistExonBinLen' => $params{'HistExonBinLen'},
	'HistIntronBinLen' => $params{'HistIntronBinLen'},
	'TrainScoreFile' => $params{'TrainScoreFile'},
	'PredScoreFile' => $params{'PredScoreFile'},
	'ConservSequence' => $params{'ConservSequence'},
	'Genes' => $params{'Genes'},
	'TrainContigs' => $params{'TrainContigs'},
	'PredContigs' => $params{'PredContigs'},
	'CrossEvalGenes' => $params{'CrossEvalGenes'},
	'Margin' => $params{'Margin'},
    };
    bless $self, $type;
    return $self;
}

sub syncInfo {
    my($self) = @_;
    if(!-e $self->{'Path'}) {
	mkdir($self->{'Path'});
    }
}

sub setGCClasses {
    my($self, $GCClasses)  = @_;
    $self->{'GCClasses'} = GenUtils::parseGCClasses($GCClasses);
    return $self->{'GCClasses'};
}

sub storeProps {
    my($self, $file) = @_;
    my $key;
    open(PROPS, ">$file") || die "Can't open $file\n";
    foreach $key (sort {$a <=> $b} keys %$self) {
	my $item;
	print PROPS "$key ";
	if($key eq 'GCClasses') {
	    my @classes = sort {$a <=> $b} keys %{$self->{$key}};
	    print PROPS join(" ", @classes), "\n";
	}
	elsif($key eq 'SignalWindows' || $key eq 'SignalModels' || $key eq 'SignalFiles' || $key eq 'IGenicGeoMeans') {
	    print PROPS join(" ", @ { $self->{$key}}), "\n";
	}
	else {
	    print PROPS "$self->{$key}\n";
	}
    }
    close(PROPS);
}

sub retrieveProps {
    my($self, $file) = @_;
    my $line = "";
    my $key;
    open(PROPS, "<$file") || die "Can't open $file\n";

    while(defined($line = <PROPS>) && $line =~ /(\S+)\s+(.+)?$/) {
	my $propName = $1; 
	my $propVal = $2; 
	if(!defined($propName) || !defined($propVal)) {
	    next;
	}
#      print "$propName, $propVal\n";
	if($propName eq 'GCClasses') {
	    $self->{'GCClasses'} = GenUtils::parseGCClasses($propVal);
	}
	elsif($propName eq 'SignalWindows' || $propName eq 'SignalModels' || $propName eq 'SignalFiles' || $propName eq 'IGenicGeoMeans') {
	    my @line = split(/\s+/, $propVal);
          $self->{$propName} = \@line;
      }
      else {
          $self->{$propName} = $propVal;
	}
    }
    close(PROPS);
}
1;
