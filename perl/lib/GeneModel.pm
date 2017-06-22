package GeneModel;
use lib "$ENV{HOME}/Lib/Perl";
use Utils;
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
	     storeProps
	     retrieveProps
	     );
sub new {
    my ($type) = shift;
    my $self = {
	'NAME' => $params{'NAME'},
	'ORG' => $params{'ORG'},
	'DOMAIN' => $params{'DOMAIN'},
	'GC_CLASSES' => $params{'GC_CLASSES'},
	'START_WINDOW' => $params{'START_WINDOW'},
	'DONOR_WINDOW' => $params{'DONOR_WINDOW'},
	'ACEPTOR_WINDOW' => $params{'ACEPTOR_WINDOW'},
	'STOP_WINDOW' => $params{'STOP_WINDOW'}
    };
    bless $self, $type;
    return $self;
}

sub storeProps {
    my($self, $file) = @_;
    open(PROPS, ">$file") || die "Can't open $file\n";
    print PROPS "$self->{'NAME'}|$self->{'ORG'}|$self->{'DOMAIN'}|$self->{'GC_CLASSES'}|$self->{'START_WINDOW'}|$self->{'DONOR_WINDOW'}|$self->{'ACEPTOR_WINDOW'}|$self->{'STOP_WINDOW'}\n";
    close(PROPS);
}

sub retrieveProps {
    my($self, $file) = @_;
    my $line = "";
    open(PROPS, "<$file") || die "Can't open $file\n";
    my $line = <PROPS>;
    defined($line) || die "Format error in $file\n";
    @line = split(/\|/, $line);
    $self->{'NAME'} = $line[0];
    $self->{'ORG'} = $line[1];    
    $self->{'DOMAIN'} = $line[2];
    $self->{'GC_CLASSES'} = $line[3];
    $self->{'START_WINDOW'} = $line[4];
    $self->{'DONOR_WINDOW'} = $line[5];
    $self->{'ACEPTOR_WINDOW'} = $line[6];
    $self->{'STOP_WINDOW'} = $line[7]; 
    close(PROPS);
}
1;
