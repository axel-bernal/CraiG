#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Utils;
use Getopt::Long;
$| = 1;
my $ds = 0;
my $ups = 0;
my $cpath = "";
my $usage = "usage : $0 -ups upstream -ds downstream -cpath path < locations\n";

&GetOptions("-ups=i" => \$ups,
            "-ds=i" => \$ds,
            "-cpath=s" => \$cpath,
	    );

(length($cpath) > 0) || die $usage;

my $cmdOpts = " -ups $ups -ds $ds ";
my $prg = "$ENV{CRAIG_HOME}/perl/bin/getSeqsFromGeneLocs.pl";
my %feats = %{&Utils::groupByPattern(\*STDIN, '^>\S+\s+(\d<|<)?([^\[,\],;]+)_\d+_\d+')};

foreach my $cid (keys %feats) {
    my $ctg = $cpath."/$cid.fa";
    if(!-e $ctg) {
        $ctg = $cpath."/$cid";
    }    
    if(!-e $ctg) {
        warn "$ctg cannot be found\n";
        return;
    }
    my $myCmdOpts = $cmdOpts." -ctg $ctg ";
    Utils::execProgram($prg, $feats{$cid}, $myCmdOpts, \*STDOUT);
}
