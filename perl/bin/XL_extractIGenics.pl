#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Utils;
use Getopt::Long;
$| = 1;
my $minSize = 1;
my $cpath = "";
my $usage = "usage : $0 -minSize minSize -cpath path < locations > iGenic_regions\n";

&GetOptions(
            "-minSize=i" => \$minSize,
            "-cpath=s" => \$cpath,
	    );

(length($cpath) > 0) || die $usage;

my $cmdOpts = " -min $minSize ";
my $prg = "$ENV{CRAIG_HOME}/perl/bin/extractIGenics.pl";
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
