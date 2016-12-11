#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Utils;
use Getopt::Long;
$| = 1;
my $usage = "usage: $0 -cpath <contig_path> < locs\n";
my $cpath = "";

&GetOptions("-cpath=s" => \$cpath, );

(length($cpath) > 0) || die $usage;

my $prg = "$ENV{CRAIG_HOME}/perl/bin/trimStopCodon.pl";
my %feats = %{&Utils::groupByPattern(\*STDIN, '^>\S+\s+(\d<|<)?([^\[,\],;]+)_\d+_\d+')};

foreach my $cid (keys %feats) {
    my $ctg = $cpath."/$cid.fa";
    print STDERR "Processing contig $ctg ..\n";
    if(!-e $ctg) {
        $ctg = $cpath."/$cid";
    }    
    if(!-e $ctg) {
        warn "$ctg cannot be found\n";
        return;
    }
    my $cmdOpts." -ctg $ctg ";
    Utils::execProgram($prg, $feats{$cid}, $cmdOpts, \*STDOUT);
}
