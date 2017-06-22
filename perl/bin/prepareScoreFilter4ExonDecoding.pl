#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(START STOP DONOR ACEPTOR);
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 scoreFile period < contigFile\n";

my $scoreFile = shift || die "Score File must be specified. ".$usage;
my $period = shift || die "period must be specified. ".$usage;
my $org = Organism->new();
my %contigs = % {$org->loadContigs(\*STDIN)};

open(SCORE, "<$scoreFile") || die "file $scoreFile not found\n";



my %signals = ("ATG", START, # start cdn -> index in 'Signals' array
               "TGA", STOP, # stop cdn
               "TAA", STOP,
               "TAG", STOP,
               "GT", DONOR, # donor signal
               "AG", ACEPTOR # aceptor signal
               );

my $line =  <SCORE>;
my $i;

while(defined($line) && $line  =~ />(\S+)\s(\d+)/) {
    my $cid = $1;
    my $clen = $2;
    my $strand = STRAND_FWD;
    
    defined($contigs{$cid}) || die "contig $cid not defined\n";
    
    my $contig = $contigs{$cid};
    $clen == $contig->{'Len'} || die "different lengths for $cid\n";    

    print ">$cid\t$clen\n";
    
    my @scores_fwd = ();
    my @scores_bwd = ();
    my $scores = [\@scores_fwd, \@scores_bwd];

    for(my $j = 0; $j <= $clen; $j++) {
        push(@{$scores->[0]}, 0);
        push(@{$scores->[1]}, 0);
    }

    while(defined($line = <SCORE>) && $line  !~ />(\S+)\s(\d+)/) {
        my $sigPos;
        my $sigScore;
        
        if($line =~ /^(\d+)\s+([-,\d,.]+)\s*$/) {
            $sigPos = $1;
            $sigScore = $2;
            $scores->[$strand]->[$sigPos] = $sigScore;                
        }
        elsif($line =~ /\/\//) {
            $strand = STRAND_COMP;
        }
    }
    

    for($strand = STRAND_FWD; $strand <= STRAND_COMP; $strand++) {
        my $seq = $contig->{'Seq'};

        if($strand == STRAND_COMP) {
            $seq = &GenUtils::revComp($seq);
        }

        my @sigArray = ();
        for(my $j = 0; $j <= $clen; $j++) {
            push(@sigArray, 0);
        }

        for($i = 1; $i <= $clen;  $i++) {
            my $signalPresent = 0;
            my $sig = uc substr($seq, $i - 1, 3);

            if(defined($signals{$sig})) {
                $signalPresent = 1;
            }
            else {
                $sig = uc substr($seq, $i - 1, 2);
                
                if(defined($signals{$sig})) {
                    $signalPresent = 1;
                }
            }

            if(!$signalPresent) {
                next;
            }

            my $k = int($i - (2*$period - 1)/2);
            if($k < 1) { $k = 1; }
            
            my $limit = int($i + (2*$period - 1)/2);
            if($limit > $clen) { $limit = $clen; }

            for( ; $k <= $limit; $k++) {
                $sigArray[$k] = 1;
                print $k,"\n";
            }
        }

        for($i = 1; $i <= $clen;  $i += $period) {
            for(my $phase = 0; $phase < $period; $phase++) {
                my $p = $i + $phase;
                if($p > $clen) {
                    next;
                }
                my $last_p = 0;
                if($p - $period > 0) {
                    $last_p = $p - $period;
                }
                $scores->[$strand]->[$p] += $scores->[$strand]->[$last_p];

                if($sigArray[$p] && $scores->[$strand]->[$p] != 0) {
                    print $p, "\t", $scores->[$strand]->[$p], "\n";
                }
            }
        }

        print "//\n";
    }
}

close(SCORE);
