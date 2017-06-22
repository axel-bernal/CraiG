#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < geneSeqFile < seqFile\n";

my $org = Organism->new();
my %contigs = % {$org->loadContigs(\*STDIN)};
my %seqs;
my $id;
my $seqId;

foreach $id (keys %contigs) {
    $id =~ /([^\@]+)\@part\:(\d+)/;
    $seqId = $1;
    my $order = $2;

    if(!defined($seqs{$seqId})) {
        @{$seqs{$seqId}} = ();
        for(my $i = 0; $i < 20; $i++) {
            push(@{$seqs{$seqId}}, "");
        }
    }

    $seqs{$seqId}->[$order] = $id;
}


foreach $seqId (keys %seqs) {
    print ">$seqId\n";
    foreach $id (@{$seqs{$seqId}}) {
        if($id eq "") {
            next;
        }
        print $contigs{$id}->{'Seq'};
    }
    print "\n";
}
