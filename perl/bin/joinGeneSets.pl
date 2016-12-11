#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : getSeqsFromLocs.pl contigFile1 locations1 contigFile2 locations2\n";

my $contigFile1 = shift || die $usage;
my $locations1 = shift || die $usage;
my $contigFile2 = shift || die $usage;
my $locations2 = shift || die $usage;

my $org1 = Organism->new();
my $org2 = Organism->new();
open(CONTIG1, "<$contigFile1") || die "file $contigFile1 not found\n";
open(LOCS1, "<$locations1") || die "file $locations1 not found\n";
open(CONTIG2, "<$contigFile2") || die "file $contigFile2 not found\n";
open(LOCS2, "<$locations2") || die "file $locations2 not found\n";
my %contigs1 = % {$org1->loadContigs(\*CONTIG1)};
close(CONTIG1);
my %genes1 = % { $org1->loadGenes(\*LOCS1)};
close(LOCS1);
my %contigs2 = % {$org2->loadContigs(\*CONTIG2)};
close(CONTIG2);
my %genes2 = % { $org2->loadGenes(\*LOCS2)};
close(LOCS2);
open(LOCS, ">$0.locs") || die;
open(CONTIGS, ">$0.contigs") || die;
my %done;
my $cid;
my $id;
foreach $id (keys %genes1) {
    $cid = $genes1{$id}->{'Exons'}->[0]->[0];
    if(defined($genes2{$id})) {
        next;
    }
    $genes1{$id}->dump(\*LOCS, \%contigs1, 1);
    print CONTIGS ">$cid\n", $contigs1{$cid}->{'Seq'}, "\n";
}

foreach $id (keys %genes2) {
    $cid = $genes2{$id}->{'Exons'}->[0]->[0];
    if(defined($genes1{$id})) {
        if($contigs1{$cid}->{'Len'} > $contigs2{$cid}->{'Len'}) {            
            $genes1{$id}->dump(\*LOCS, \%contigs1, 1);
            print CONTIGS ">$cid\n", $contigs1{$cid}->{'Seq'}, "\n";
            next;
        }
    }
    $genes2{$id}->dump(\*LOCS, \%contigs2, 1);
    print CONTIGS ">$cid\n", $contigs2{$cid}->{'Seq'}, "\n";            
}

close(CONTIGS);
close(LOCS);
