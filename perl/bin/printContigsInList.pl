#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use GenUtils;
use Organism;
use Contig;
use Gene;
use strict;
use Getopt::Long;
$| = 1;

# -*- perl -*-

my $usage = "usage: $0 contigs < list of contig ids > contigs in list\n";
 
my $org = Organism->new();
my $contigFile = shift || die $usage;
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
my $id;
my $line = "";

while(defined($line = <STDIN>))  {
    $line =~ /^(\S+)\s+/;
    $id = $1;
    chomp($line);
    defined($contigs{$id}) || next;
    print ">$line\n$contigs{$id}->{'Seq'}\n";
}
