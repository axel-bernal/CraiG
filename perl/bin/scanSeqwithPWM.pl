#!/usr/bin/perl -w
#Scan for matches of position weight matrix over a set of sequences
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use GenUtils;
use PWM;
use strict;
use Getopt::Long;
use File::Basename;

my $line = "";
my $prof = "";
my $filter = "";
my $cutoff = 10;
my $strand = 0;
my $pwm_mats;
my $pwm_filters;
my $filt_spacers;
my $spacers;
my @positions = ();

my $i;
my $j;

my $usage = "$0 --pwm profile -filter filter --cutoff cutoff --strand strand < sequence\nProfile contains matrix with the weights for each position and character found in sequence";

&GetOptions(
            "pwm=s" => \$prof,
	    "filter=s" => \$filter,
            "cutoff=f" => \$cutoff,
	    "strand=i" => \$strand,
	    );

length($prof) != 0 || die $usage;
open(PROF,"$prof") || die "Cannot find/open profile $prof\n";

# reading the matrix information
my $pwm = new PWM();
$pwm->readPWM(\*PROF,'\s+');
$pwm_mats = $pwm->{'Matrices'};
$spacers = $pwm->{'DisToNextMat'};
#foreach $i (keys % { $pwm_mats->[1]}) {
#   print join(" ",@ {$pwm_mats->[1]->{$i}}),"\n";
#}

#reading filter information
if(length($filter) > 0) {
    open(FILTER, "<$filter") || die "Can't open $filter\n";
    $pwm->setFilter(\*FILTER,'\s+');
}
else {
    $pwm->setFilter();
}
$pwm_filters = $pwm->{'Filters'};

# scanning it through each file
$line = <>;
my ($fp, $tp, $fn) = (0,0,0,0);
while(defined($line) && $line =~ /^>/) {
    my $header = $line;
    my @seq = ();
    my $seq;
    my $i = 0;
    while(defined($line = <>) && $line !~ /^>/) {
	push(@seq, $line);
    }
    map {chomp} @seq;
    $seq = join("", @seq);
    $seq = uc $seq;
    print $header;
    # hack added must remov later
    $header =~ /\S+\s+(.+)$/;

    my (undef, undef, $exons, undef) = @{&GenUtils::readContigLocation($1)};
    if(!$strand || $strand == 1) {
	@positions = @ {&PWMUtils::signalPosFromPWM($pwm, $seq, $cutoff)};
	print "#forward strand :" , scalar(@$pwm_mats)," pwms\n";
	my $e = 1;
	my $thisTP = 0;
	for($j = 0; $j < scalar(@ positions);$j++) {
	    if(scalar(@$exons) > $e) {
		if($positions[$j]->[0] == $exons->[$e]->[1] - 15) {$thisTP++; $e++;}
	    }
	    print $positions[$j]->[0]," ",$j;
	    for($i = 0;$i < scalar(@$pwm_mats);$i++) {
		print " ", $positions[$j]->[$i+1], " ", substr($seq,$positions[$j]->[$i],scalar(@{$pwm_mats->[0]->{chars}}));
	    }
	    print "\n";
	}
	$fn += scalar(@$exons) - 1 - $thisTP;
	$fp += scalar(@positions) - $thisTP;
	$tp += $thisTP
    }
    if(!$strand || $strand == -1)    {
	my $rev_seq = &GenUtils::revComp($seq);
	@positions = @ {&PWMUtils::signalPosFromPWM($pwm, $rev_seq, $cutoff)};
	print "#backward strand :", scalar(@$pwm_mats)," pwms\n";
	for($j = 0; $j < scalar(@ positions);$j++) {
	    for($i = 0;$i < scalar(@ {$positions[$j]});$i++) {
		print " $positions[$j]->[$i]";
	    }
	    for($i = 0;$i < scalar(@$pwm_mats);$i++) {
		print " ",substr($rev_seq,$positions[$j]->[$i],scalar(@{$pwm_mats->[0]->{chars}}));
	    }
	    print "\n";
	}
    }
}
print "$tp $fp $fn ", $tp/($tp+$fp)," ", $tp/($tp + $fn),"\n";
