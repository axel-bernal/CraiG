package Histogram;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
use Utils;
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
	     smooth
	     smoothedValue
	     value
	     );

use strict;

sub new {
    my $type  = shift;
    my %params = @_;
    my $self = { 
	'OriginValues' => $params{'OriginValues'},
	'SmoothValues' => undef,
	'SortedKeys' => undef,
	'Min' => undef,
	'Max' => undef,
    };
    @ { $self->{'SortedKeys'}} = sort {$a <=> $b} (keys %{$self->{'OriginValues'}});
    if(scalar(@ { $self->{'SortedKeys'}}) == 0) {
	$self->{'Min'} = $self->{'Max'} = -1;
    }
    else {
	$self->{'Min'} = $self->{'SortedKeys'}->[0];
	$self->{'Max'} = $self->{'SortedKeys'}->[scalar(@ { $self->{'SortedKeys'}}) - 1];
    }
    bless $self, $type;
    return $self;
}

sub smooth {
    my ($self) = @_;
    my $tmpFile = "$ENV{CRAIG_HOME}/tmp/$$.hist";
    if(scalar(@{$self->{'SortedKeys'}}) <= 2) {
	for(my $i = $self->{'Min'}; $i <= $self->{'Max'}; $i++) {
	    $self->{'SmoothValues'}->{$i} = $self->{'OriginValues'}->{$self->{'Min'}} + ($self->{'OriginValues'}->{$self->{'Max'}} - $self->{'OriginValues'}->{$self->{'Min'}})*$i;
	}
	return $self->{'SmoothValues'};
    }
    &Utils::hashToFile($self->{'OriginValues'},"$tmpFile.in"," ");
    open(FILE, ">$tmpFile") || die "Can't open $tmpFile\n";
    print FILE "set output \'$tmpFile.out\'\n";
    print FILE "set term table\n";
    print FILE "set samples ",$self->{'Max'} - $self->{'Min'} + 1,"\n";
    print FILE "plot \'$tmpFile.in\' s c\n";
    close(FILE);
    `gnuplot $tmpFile`;
    unlink($tmpFile);
    unlink("$tmpFile.in");
    $self->{'SmoothValues'} = Utils::fileToHash("$tmpFile.out", " ", 1);
    foreach my $val (sort {$a <=> $b} keys %{$self->{'SmoothValues'}}) {
	if($self->{'SmoothValues'}->{$val} eq "" || $self->{'SmoothValues'}->{$val} < 1) {
	    $self->{'SmoothValues'}->{$val} = 1;
	}
    }
    # minimum length is 3 = codon
    undef($self->{'SmoothValues'}->{1});
    undef($self->{'SmoothValues'}->{2});
    unlink("$tmpFile.out");
    return $self->{'SmoothValues'};
}

sub smoothedValue {
    my($self, $pos) = @_;
    defined($self->{'SmoothValues'}->{$pos}) || die "histogram value not found at $pos\n";
    return $self->{'SmoothValues'}->{$pos};
}

sub value {
    my($self, $pos) = @_;
    defined($self->{'OriginValues'}->{$pos}) || die "histogram value not found at $pos\n";
    return $self->{'OriginValues'}->{$pos};
}
1;
