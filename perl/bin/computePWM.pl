#!/usr/bin/perl -w
#Build a PWM from a set of sewuences
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use GenUtils;
use PWM;
use strict;
use Getopt::Long;
use File::Basename;

my $pwm_mats;
my $pwm_filters;
my $filt_spacers;
my $spacers;
my @positions = ();

my $i;
my $j;

my $usage = "$0 <File with background frequency> < sequence > PWM \nPWM contains matrix with the weights for each position and character found in sequence";
my $bkgFile = shift @ARGV || die $usage;
open(BKG, "$bkgFile") || die;
# reading the matrix information
my $pwm = new PWM();

sub readFasta {
    my $file = shift @_;
    my $line = "";
    my @seqs = ();
    $line = <$file>;
    while(defined($line) && $line =~ /^>/) {
	my $header = $line;
	my @seq = ();
	my $seq;
	my $i = 0;
	while(defined($line = <$file>) && $line !~ /^>/) {
	    push(@seq, $line);
	}
	map {chomp} @seq;
	$seq = join("", @seq);
	$seq = uc $seq;
	$seq =~ tr/NUMRWSYKBDHV/AAAAAAAAAAAA/;
	push(@seqs, $seq);
    }
    return \@seqs;
}

my @seqs = @ {&readFasta(\*STDIN)};
my @bkgSeqs = @ {&readFasta(\*BKG)};


$pwm_mats = $pwm->{'Matrices'};
$spacers = $pwm->{'DisToNextMat'};
#foreach $i (keys % { $pwm_mats->[1]}) {
#    print join(" ",@ {$pwm_mats->[1]->{$i}}),"\n";
#}
$pwm->computePWM(\@seqs, \@bkgSeqs);
$pwm->printPWMFmt(\*STDOUT);
