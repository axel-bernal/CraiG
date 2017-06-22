package PWMUtils;
require Exporter;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
use Utils;
@ISA = ( Exporter );
@EXPORT= qw (
	     scan4PWMMatches
	     signalPosFromPWM
	     );

sub scan4PWMMatches {
    my($seq, $pwmObj, $i) = @_;
    my ($pwm, $filter) = ($pwmObj->{'Matrices'}->[$i], $pwmObj->{'Filters'}->[$i]);
    my @results = ();
    my $pos = 0;
    for($i = 0; $i < length($seq);$i++) {
        $results[$i] = 0;
    }
    while($pos + $pwm->{"chars"}->[scalar(@ {$pwm->{"chars"}}) - 1] <= length($seq))    {        
	for($i = 0; $i < scalar(@ {$pwm->{"chars"}}); $i++) {
            $char = substr($seq,$pos + $pwm->{"chars"}->[$i] - 1,1); 
            if(!defined($pwm->{$char})) {
                next;
            }
            if($filter->{$char}->[$i] == 0) {
                $results[$pos + ($pwm->{"chars"})->[0] - 1] = 0;
                last;
            }
            else {
                $results[$pos + ($pwm->{"chars"})->[0] - 1] += ($pwm->{$char})->[$i];
	    }
	}
#	print $pos, " | " ,$pwm->{"chars"}->[scalar(@ {$pwm->{"chars"}}) - 1]," | ", $results[$pos + ($pwm->{"chars"})->[0] - 1],"\n";
	$pos++;
    }    
    return \@results;
}

# Compute the positions of occurrances of signals best described by the PWM in sequence seq. Only reporting those above cutoff
sub signalPosFromPWM {
    my ($pwm, $seq, $cutoff) = @_;
    my($pwm_mats, $pwm_filters, $spacers) =  ($pwm->{'Matrices'}, $pwm->{'Filters'}, $pwm->{'DisToNextMat'});
    my @scores = ();
    my @mat_scores = ();
    my @bt_matrix = ();
    my(@results) = ();
    my $acc_offset = 0;
    my $max_score = 0;
    my $acc_min_offsets = 0;
    my ($i, $j, $st) = (0,0,0);
    for($i = 0; $i < scalar(@{$pwm_mats});$i++) {
        @ { $scores[$i]} = @ {&PWMUtils::scan4PWMMatches($seq, $pwm, $i)};
#        print " ", $i , " " ,scalar(@{$pwm_mats}), " ", join(" ", @ { $scores[$i]}), "\n";
    }
    # initializing each matrix
    for($i = 0; $i < scalar(@{$pwm_mats}); $i++) {
	for($j = 0; $j < scalar(@ { $scores[$i]}); $j++) {
	    if($i == 0) {
		$mat_scores[$i]->[$j] = $scores[$i]->[$j];
		$bt_matrix[$i]->[$j] = 0;
	    }
	    else {
		$mat_scores[$i]->[$j] = 0;
		$bt_matrix[$i]->[$j] = -1;              
	    }
	}
#       print "AAAA ",join("\n",@{$scores[$i]}),"\n";
    }
    
    # iteration through each of the matrix scores
    for($st = 1; $st < scalar(@{$pwm_mats}); $st++) {
	my $len_signal = (($pwm_mats->[$st - 1])->{"chars"})->[scalar(@{($pwm_mats->[$st - 1])->{"chars"}}) - 1];
	my $min_offset = $spacers->[$st - 1]->[0] + $len_signal - 1;
	my $max_offset = $spacers->[$st - 1]->[1] + $len_signal - 1;
	$acc_offset += $max_offset + $acc_min_offsets;
	$acc_min_offsets += $min_offset;
#       print "STATE = $st $len_signal $min_offset $max_offset $acc_offset\n";
	for($j = $acc_offset; $j < scalar(@ { $scores[$st]}); $j++) {
	    $max_score = 0;
#           print "comming = $j\n";
	    for($i = $j - $max_offset; $i <= $j - $min_offset; $i++) {
		if($bt_matrix[$st - 1]->[$i] < 0) {
		    next;
		}
#               print "old = $i\n";
		if($mat_scores[$st - 1]->[$i] + $scores[$st]->[$j] > $max_score) {
                    $mat_scores[$st]->[$j] = $mat_scores[$st - 1]->[$i] + $scores[$st]->[$j];
		    $bt_matrix[$st]->[$j] = $i;
		    $max_score = $mat_scores[$st]->[$j];
		}
	    }
#           print "maximum $j $i = ", $bt_matrix[$st]->[$j]," score = $mat_scores[$st]->[$j]\n";
        }
    }
    # get only those above cutoff value
    for($i = 0; $i < scalar(@ { $mat_scores[scalar(@{$pwm_mats}) - 1]}); $i++) {
        my @bt_path = ();
	my @ind_scores = ();
        if(($mat_scores[scalar(@{$pwm_mats}) - 1])->[$i] < $cutoff) {
            next;
        }
        $bt_path[scalar(@{$pwm_mats}) - 1] = $i;
        $ind_scores[scalar(@{$pwm_mats})] = ($mat_scores[scalar(@{$pwm_mats}) - 1])->[$i];
        # backtracking
        for($j = scalar(@{$pwm_mats}) - 1; $j > 0 ; $j--) {
            $ind_scores[$j] = $scores[$j]->[$bt_path[$j]];
            $bt_path[$j - 1] = $bt_matrix[$j]->[$bt_path[$j]];
        }
        $ind_scores[0] = $scores[0]->[$bt_path[0]];
        push(@bt_path, @ind_scores);
        push(@results, \@bt_path);
    }
    return \@results;
}
