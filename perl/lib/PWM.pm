package PWM;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
use PWMUtils;
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
	     readPWM
	     computePWM
	     );

sub new {
    my $type = shift;
    my %params = @_;
    my $self = {
	'Matrices' => $params{'Matrices'},
	'DisToNextMat' => $params{'DisToNextMat'},
	'Filters' => $params{'Filters'},
    };
    bless $self, $type;
    return $self;
}

# Read a PWM structure. Can be composed of more than one entry separated by Spacers which are variable distance separators. Here is example :                  
#chars 1 2 3 4 5 6 7 8 9 10                                                  
#A 14 11 26 0 28 2.5 0 0.5 0 3                                                
#C 3 1 1 36 1 0.5 0 24.5 33 23                                                 
#G 16 24 10 0 0 0 37 0 0 0                                                     
#T 4 1 0 1 7 34 0 12 4 10                                                      
#[0..14]                                                                       
#chars 11 12 13 14 15 16 17 18 19 20                                           
#A 6 2 11.5 0 27 4 0 0.5 1 2                                                   
#C 2 0 0.5 36 2 0 0 9.5 24 15                                                  
#G 23.5 34 25 0 2 1 37 0 0 3                                                   
#T 5.5 1 0 1 5 32 0 27 12 16                                                   
#//                                 
sub readPWM {
    my ($self, $file, $col_sep) = @_;
    my @spacers = ();
    my @pwms = ();
    my $count = 0;
    my($line, $row_key);
    %{$pwms[$count]} = ();
    while(defined($line = <$file>) && $line !~ /\/\//) {
	if($line =~ /\[(\d+)\.\.(\d+)\]/) { # spacer                          
	    push(@spacers, [$1, $2]);
	    %{$pwms[++$count]} = ();
	    next;
	}
	my(@line) = split($col_sep,$line);
	$row_key = $line[0];
	shift @line;
	push(@ {$pwms[$count]->{$row_key}},@line);
    }
    $self->{'Matrices'} = \@pwms;
    $self->{'DisToNextMat'} = \@spacers;
    return $self;
}

sub setFilter {
    my $self = shift @_;
    my ($file, $sepCol) = @_;
    my $distances;
    if(defined($file)) {
	my  $tmp = new PWM();
	$tmp->readPWM($file, $sepCol);
	$self->{'Filters'} = $tmp->{'Matrices'}; 
	$distances = $tmp->{'DisToNextMat'};
    }
    else {
	my $i;
	my $pwm_filters;
	for($i = 0; $i < scalar(@{$self->{'Matrices'}});$i++) {
	    %{$pwm_filters->[$i]} = ();
	    my $row_key;
	    foreach $row_key (keys %{$self->{'Matrices'}->[$i]}) {
		my $j;
		for($j = 0; $j < scalar(@ {$self->{'Matrices'}->[$i]->{$row_key}});$j++) {
		    if($row_key eq "chars") {
			$pwm_filters->[$i]->{$row_key}->[$j] = $self->{'Matrices'}->[$i]->{$row_key}->[$j];
		    }
		    else {
			$pwm_filters->[$i]->{$row_key}->[$j] = 1;
		    }
		}
	    }
	}
	$self->{'Filters'} = $pwm_filters;
	$distances = $self->{'DisToNextMat'};
    }
    (!defined($distances) || scalar(@{$distances}) == scalar(@{$self->{'DisToNextMat'}})) || die "Incompatible filter\n";
    (scalar(@{$self->{'Matrices'}}) == scalar(@{$self->{'Filters'}})) || die "Incompatible filter\n";    
}

sub computePWM {
    my($self, $seqs, $bkgSeqs) = @_;
    my %expFreqs = ();
    my %posFreqs = ();
    my $total = 0;
    my @posTotal = ();
    my $i;
    my @pwms = ();
    my $count = 0;
    %{$pwms[$count]} = ();
    my $c;
    (scalar(@$seqs) != 0) || die;
    my $len = length($seqs->[0]);
    my $seq;
    foreach $seq (@$seqs) {
	for($i = 0; $i < $len; $i++) {
	    $c = substr($seq, $i, 1);
	    $posFreqs{$c}->[$i]++;
	    $posTotal[$i]++;
	}
   }
    foreach $seq (@$bkgSeqs) {
	my $len = length($seq);
	for($i = 0; $i < $len; $i++) {
	    $c = substr($seq, $i, 1);
	    $expFreqs{$c}++;
	    $total++;
	}
    }

    for($i = 0; $i < $len; $i++) {
	$pwms[$count]->{'chars'}->[$i] = $i + 1;
    }
    
    foreach $c (keys %posFreqs) {
	for($i = 0; $i < $len; $i++) {
	    if(!defined($posFreqs{$c}->[$i])) { 
		$pwms[$count]->{$c}->[$i] = log((0.0000000000000001*$total)/($expFreqs{$c}*$posTotal[$i]))/ log(2);
	    }
	    else { 
		$pwms[$count]->{$c}->[$i] = (log($posFreqs{$c}->[$i]/$posTotal[$i]) - log($expFreqs{$c}/$total))/ log(2);
	    }
	}
    }
    $self->{'Matrices'} = \@pwms;
    return $self;
}

sub printPWMFmt {
    my($self, $file) = @_;
    return unless defined($self->{'Matrices'});  
    my ($i, $j);
    my $row_key;
# pwm values
    for($i = 0; $i < scalar(@{$self->{'Matrices'}});$i++) {
	if($i > 0) {
	    print $file "[", join("..", @{$self->{'DisToNextMat'}->[$i]}),"]\n";
	}
	print $file "chars ", join(" ", @ {$self->{'Matrices'}->[$i]->{chars}}), "\n";
	foreach $row_key (keys %{$self->{'Matrices'}->[$i]}) {
	    if($row_key eq "chars") {
		next;
	    }
	    print $file $row_key, " ";
	    for($j = 0; $j < scalar(@ {$self->{'Matrices'}->[$i]->{$row_key}}); $j++) {
		print sprintf "%3.2f",$self->{'Matrices'}->[$i]->{$row_key}->[$j];
		print " ";
	    }
	    print "\n";
#	    print $file $row_key, " ", join(" ", @ {$self->{'Matrices'}->[$i]->{$row_key}}), "\n";
	}
    }	
    print $file "//\n";
    return unless defined($self->{'Filters'});  
# filters
    for($i = 0; $i < scalar(@{$self->{'Filters'}});$i++) {
	if($i > 0) {
	    print $file "[", join("..", @{$self->{'DisToNetMat'}->[$i]}),"]\n";
	}
	print $file "chars ", join(" ", @ {$self->{'Filters'}->[$i]->{chars}}), "\n";
	foreach $row_key (keys %{$self->{'Filters'}->[$i]}) {
	    if($row_key eq "chars") {
		next;
	    }
	    print $file $row_key, " ", join(" ", @ {$self->{'Filters'}->[$i]->{$row_key}}), "\n";
	}
    }	
    print $file "//\n";
}
