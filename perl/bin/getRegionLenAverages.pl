#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
use Organism;
use Contig;
use Getopt::Long;

my $usage = "usage : $0 contigFile < locations\n\tOutputs length averages for introns, exons, igenics\n";

my $contigFile = shift || die $usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my %genes = % { $org->loadGenes(\*STDIN)};

my @initials = ();
my @internals = ();
my @terminals = ();
my @singles = ();
my @introns = ();
my @igenics = ();

for my $c (keys %contigs) {
    my $genes = $contigs{$c}->{'Features'}->{'TRAIN'};

    foreach my $id (keys %$genes) {
        my $numExons = scalar(@{$genes->{$id}->{'Exons'}});
        for(my $i = 0; $i < $numExons; $i++) {
            my $ce = $genes->{$id}->{'Exons'}->[$i];
            if($i > 0) {
                my $le = $genes->{$id}->{'Exons'}->[$i - 1];
                if($genes->{$id}->{'Strand'} == STRAND_FWD) {
                    push(@introns, [$ce->[0], $le->[2] + 1, $ce->[1] - 1]);
                }
                else {
                    push(@introns, [$ce->[0], $le->[2] - 1, $ce->[1] + 1]); 
                }
            }

            if($numExons == 1) {  push(@singles, $ce); }
            else {
                if($i == 0) {  push(@initials, $ce); }
                elsif($i == $numExons - 1)  {  push(@terminals, $ce); }
                else  {  push(@internals, $ce); }
            }                                           
        }
    }
    
    # dealing with igenics
    my $array = &GenUtils::markIGenics($contigs{$c}, 0, 1, 'TRAIN');
    my $interBeg = -1;
    my $interEnd = -1;

    for(my $i = 1; $i <= $contigs{$c}->{'Len'}; $i++) {
        if($array->[$i] == 0 && $interBeg == -1) {
            $interBeg = $i;
        }
        if($i > 1) {
            if(($array->[$i] == 1 && $array->[$i - 1] == 0) || $i == $contigs{$c}->{'Len'}) {
                $interEnd = $i - 1;
            }
        }

        if($interBeg != -1 && $interEnd != -1) {
            push(@igenics, [$c, $interBeg, $interEnd]);
            $interBeg = -1;
            $interEnd = -1;
        }        
    }
}

print $#initials + 1, " initial exons\n";
printLengthAverages(\@initials);

print $#internals + 1, " internal exons\n";
printLengthAverages(\@internals);

print $#terminals + 1, " terminal exons\n";
printLengthAverages(\@terminals);

print $#singles + 1, " single exons\n";
printLengthAverages(\@singles);

print $#introns + 1, " introns\n";
printLengthAverages(\@introns);

print $#igenics + 1, " intergenics\n";
printLengthAverages(\@igenics);

sub printLengthAverages {
    my($array) = @_;
    my $total = 0;

    for(my $i = 0; $i < $#$array; $i++) {
	my $e = $array->[$i];
        $total += abs($e->[2] - $e->[1]) + 1;
    }
    
    print "average = ", $total/($#$array + 1), "\n";
}
