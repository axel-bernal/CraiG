#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Utils;
use Getopt::Long;
$| = 1;
my $outC = "";
my $trimVal = 2000;
my $cpath = "";
my $locs = 0;
my $usage = "usage : $0 [-l] -d trimVal -outC out_contig -cpath path < locations. If -l is specified, outContigs will have the new contigs as locations.\n";

&GetOptions("-outC=s" => \$outC,
            "-d=i" => \$trimVal,
            "-cpath=s" => \$cpath,
            "-l" => \$locs
	    );

(length($cpath) > 0 && length($outC) > 0) || die $usage;

my $cmdOpts = ($locs ? " -l " : " ");
$cmdOpts .= "-d $trimVal ";
my $prg = "$ENV{CRAIG_HOME}/perl/bin/trimGeneSeqs.pl";
my %feats = %{&Utils::groupByPattern(\*STDIN, '^>\S+\s+(\d<|<)?([^\[,\],;]+)_\d+_\d+')};
unlink($outC);

foreach my $cid (keys %feats) {
    my $ctg = $cpath."/$cid.fa";
    if(!-e $ctg) {
        $ctg = $cpath."/$cid";
    }
    if(!-e $ctg) {
        warn "$ctg cannot be found\n";
        return;
    }
    my $outCtg = "$outC.$$";

    my $myCmdOpts = $cmdOpts." -ctg $ctg -outC $outCtg ";
    Utils::execProgram($prg, $feats{$cid}, $myCmdOpts, \*STDOUT);
    `cat $outCtg >> $outC`;
    unlink($outCtg);
}
