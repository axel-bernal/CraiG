#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $usage = "usage : $0 contigFile < fasta\n";

$contigFile = shift || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my $line;
while(defined($line = <STDIN>)) {
    if($line =~ /^>(\S+)/) {
        chomp $line;
        my $id = $1;
        if(defined($contigs{$id})) {
            print $line, "\t", $contigs{$id}->{'Len'}, "\n";
        }
        else { die "bad cid  $id\n"};
    }
    else { print $line;}
}
