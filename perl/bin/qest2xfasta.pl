#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use GenUtils;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $usage = "usage : $0  marked_coords qests < subcontig_info\n";

my $markFile = shift || die $usage;
my $estsFile = shift || die $usage;
my $org = Organism->new();
my $org2 = Organism->new();
my $id;
my $i;
my $mark;

open(MARKS, "<$markFile") || die "file $estsFile not found\n";
open(ESTS, "<$estsFile") || die "file $estsFile not found\n";
my %marks = %{&GenUtils::loadExtendedFasta(\*MARKS)};
my %ests = %{&GenUtils::loadExtendedFasta(\*ESTS)};
my %estSeqs = ();

foreach $id (keys %ests) {
    my $seq = "";
    for($i = 0; $i < scalar(@{$ests{$id}}); $i++) {
        if($ests{$id}->[$i]->[1]) {
            $seq = $seq.($ests{$id}->[$i]->[0] x $ests{$id}->[$i]->[1]);
        }
        else {
            print STDERR "0 length: $id $ests{$id}->[$i]->[0] $ests{$id}->[$i]->[1] $i\n";
        }
    }
    $estSeqs{$id} = $seq;
}

my %subcontigs = % {$org2->loadContigs(\*STDIN)};
my %estxSeqs = ();
foreach $id (keys %subcontigs) {
    $id =~ /^(\S+)\.(\d+)$/;
#    $id =~ /^(\S+)$/;
    my $len = $subcontigs{$id}->{'Len'};
    my $estLen = length($estSeqs{$1."+"});
    my $_seqPlus = substr($estSeqs{$1."+"}, $2 - 1, $len);
    my $_seqMinus = substr($estSeqs{$1."-"}, $estLen - ($2 + $len - 1), $len);
    my $seqPlus = GenUtils::fasta2xfasta($_seqPlus);
    my $seqMinus = GenUtils::fasta2xfasta($_seqMinus);
    $estxSeqs{$id} = [$seqPlus, $seqMinus];
    print ">$id\n$estxSeqs{$id}->[0]\/\/\n$estxSeqs{$id}->[1]";
}
exit;

foreach $id (keys %estxSeqs) {
    my $mark = $marks{$id};
    if($mark->[0]->[0] == 0) {
        $estxSeqs{$id}->[0] = "- x $mark->[0]->[1]\n".$estxSeqs{$id}->[0];
        $estxSeqs{$id}->[1] = $estxSeqs{$id}->[1]."- x $mark->[0]->[1]"."\n";
    }
    if(scalar(@{$mark}) > 1 && $mark->[-1]->[0] == 0) {
        $estxSeqs{$id}->[0] = $estxSeqs{$id}->[0]."- x $mark->[-1]->[1]"."\n";
        $estxSeqs{$id}->[1] = "- x $mark->[-1]->[1]\n".$estxSeqs{$id}->[1];
    }

    print ">$id\n$estxSeqs{$id}->[0]\/\/\n$estxSeqs{$id}->[1]";
    
}

