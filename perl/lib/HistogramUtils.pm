package HistogramUtils;
use lib "$ENV{CRAIG_HOME}/Lib/Perl";
use Utils;
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
	     initGenomicRegLengths
	     extractGenomicRegLengths
             getLengthControlPoints
             printLengthControlPoints
             createLengthModelParamsFile
             computeBinningSizes
	     );

use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

sub initGenomicRegLengths {
    my ($stateLens, $numCPs, $cpExons, $cpIntrons, $has_lintron) = @_;
    $has_lintron = 0 if !defined($has_lintron);
    my $intron_st = "INTRON";
    $intron_st = "LONG_".$intron_st if $has_lintron;
    @{$stateLens->{"INIT_EXON"}} = ();
    @{$stateLens->{"INTERNAL_EXON"}} = ();
    @{$stateLens->{"LAST_EXON"}} = ();
    @{$stateLens->{"SINGLE_EXON"}} = ();
    @{$stateLens->{"INIT_UTR"}} = ();
    @{$stateLens->{"INTERNAL_UTR"}} = ();
    @{$stateLens->{"LAST_UTR"}} = ();
    @{$stateLens->{"SINGLE_UTR"}} = ();

    @{$stateLens->{$intron_st}} = ();
    @{$stateLens->{"UTR_INTRON"}} = ();

    $numCPs->{"INIT_EXON"} = $cpExons;
    $numCPs->{"INTERNAL_EXON"} = $cpExons;
    $numCPs->{"LAST_EXON"} = $cpExons;
    $numCPs->{"SINGLE_EXON"} = $cpExons;
    $numCPs->{"INIT_UTR"} = $cpExons;
    $numCPs->{"INTERNAL_UTR"} = $cpExons;
    $numCPs->{"LAST_UTR"} = $cpExons;
    $numCPs->{"SINGLE_UTR"} = $cpExons;

    $numCPs->{$intron_st} = $cpIntrons;
    $numCPs->{"UTR_INTRON"} = $cpIntrons;
}

sub extractGenomicRegLengths {
    my($stateLens, $exons, $strand, $has_lintron) = @_;
    $has_lintron = 0 if !defined($has_lintron);
    my $intron_st = "INTRON";
    $intron_st = "LONG_".$intron_st if $has_lintron;
    my $numExons = scalar(@$exons);
    for(my $i = 0; $i < $numExons; $i++) {
	my $ce = $exons->[$i];
	if($i > 0) {
	    my $le = $exons->[$i - 1];
	    if($strand == STRAND_FWD) {
		my $intron = [$ce->[0], $le->[2] + 1, $ce->[1] - 1];
		push(@{$stateLens->{$intron_st}}, $intron);
		push(@{$stateLens->{"UTR_INTRON"}}, $intron);
	    }
	    else {
		my $intron = [$ce->[0], $le->[2] - 1, $ce->[1] + 1];
		push(@{$stateLens->{$intron_st}}, $intron);
		push(@{$stateLens->{"UTR_INTRON"}}, $intron);
	    }
	}
	
	if($numExons == 1) {  
	    push(@{$stateLens->{"SINGLE_EXON"}}, $ce);
	    push(@{$stateLens->{"SINGLE_UTR"}}, $ce);
	}
	else {
	    if($i == 0) {
		push(@{$stateLens->{"INIT_EXON"}}, $ce);
		push(@{$stateLens->{"INIT_UTR"}}, $ce);
	    }
	    elsif($i == $numExons - 1)  {
		push(@{$stateLens->{"LAST_EXON"}}, $ce);
		push(@{$stateLens->{"LAST_UTR"}}, $ce);
	    }
	    else  {  
		push(@{$stateLens->{"INTERNAL_EXON"}}, $ce);
		push(@{$stateLens->{"INTERNAL_UTR"}}, $ce);
	    }
	}                                           
    }
}

sub getLengthControlPoints {
    my($st_lens, $num_cps) = @_;
    my $offset = 0;
    my @ordered = sort {abs($a->[2] - $a->[1]) <=> abs($b->[2] - $b->[1]); } @$st_lens;
    my $e;
    my @cps = ();
    for(my $i = 0; $i < $num_cps; $i++) {
        $e = $ordered[int($offset)];
	push(@cps, [int($offset), abs($e->[2] - $e->[1]) + 1]);
        $offset += 1.0*scalar(@ordered)/$num_cps;
    }
    $e =  $ordered[-1];
    push(@cps, [$#ordered, abs($e->[2] - $e->[1]) + 1]);
    return \@cps;
}

sub printLengthControlPoints {
    my($file, $st_lens, $cps, $keyw) = @_;
    print $#$st_lens + 1, " $keyw states\n";
    for(my $i = 0; $i <= $#$cps; $i++) {
        print $file "[", $cps->[$i]->[0],"] -> ", $cps->[$i]->[1], "\n";
    }
}

sub createLengthModelParamsFile  {
    my($cps, $tparams, $oparams) = @_;
    open(IN, $tparams) || die "Can't open $tparams\n";
    open(OUT, ">$oparams") || die "Can't open $oparams\n";
    my $line;

    while(defined($line = <IN>)) {
	if($line =~ /^Lengths(\s+)(\d?)(\S+)\_(\S)(\s+)\[([^\]]+)\]/) {
	    my $ranges = $6;
	    my $state = $3;
	    $line = "Lengths$1$2$3_$4$5\[";

	    if(defined($cps->{$state})) {
		while($ranges =~ /^(\s*\d+)\:(\d+)\-\>([\-,\d]+)(.*)$/) {
		    if($3 != -1) {
			$line .= "$1"."\:"."$2"."->".int(1.05*$cps->{$state}->[-1]->[1]);
		    }
		    else {
			$line .= "$1"."\:"."$2"."->$3";
		    }
		    $ranges = $4;
		}
		$line .= "\]\n";
	    }
	    else {
		$line = $line."$ranges\]\n";
	    }
	}
	print OUT $line;
    }
    
    close(IN);
    close(OUT);
}


sub computeBinningSizes {
    my($subs, $cps) = @_;
    foreach my $keyw (keys %$cps) {
	die "no $keyw states found\n" if !scalar(@{$cps->{$keyw}});
	my @lenvals = ();

	for(my $i = 0; $i < scalar(@{$cps->{$keyw}}); $i++) {
	    push(@lenvals, $cps->{$keyw}->[$i]->[1]);
	}

	my $best_step = 2;
	my $best_mean_ratio = 1;
	my $best_numbins = 64;
	my $best_stdev = 1000000;
#	print "state $keyw\n";
	for(my $step = 6; $step < 13; $step++) {
	    my @ratios = ();
	    my $num_bins = 0;
	    my $i = $step - 1;
	    my $next_i = $i + $step;
	    my $last_val = 0;
	    my $mean_ratio = 0;
#	    print "step $step\n";
	    while($next_i < 64) {
		push(@ratios, ($lenvals[$next_i] - $lenvals[$i])/($lenvals[$i] - $last_val));
		$mean_ratio += $ratios[-1];
#		print "ratio $lenvals[$next_i], $lenvals[$i], $last_val, $next_i $ratios[-1]\n";
		$last_val = $lenvals[$i];
		$i = $next_i;
		$next_i = $i + $step;
		$num_bins++;
	    }

	    $mean_ratio /= $num_bins;

	    my $ratio_stdev = 0;
	    for($i = 0; $i < $num_bins; $i++) {
		$ratio_stdev += ($ratios[$i] - $mean_ratio)**2;
	    }
	    $ratio_stdev /= $num_bins;

#	    print "step $step $mean_ratio $ratio_stdev $num_bins\n";

	    if($ratio_stdev < $best_stdev) {
		$best_stdev = $ratio_stdev;
		$best_numbins = $num_bins;
		$best_mean_ratio = $mean_ratio;
		$best_step = $step;
	    }
	}

	$subs->{$keyw."_NUMBINS"} = $best_numbins + 1;
	$subs->{$keyw."_LOGBASE"} = sprintf("%.1f", $best_mean_ratio);
	$subs->{$keyw."_BINSIZE"} = $lenvals[$best_step];
    }
}
