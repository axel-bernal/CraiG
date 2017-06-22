#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Utils;
use Getopt::Long;
$| = 1;
my $dirC = "";
my $minSize = 0;
my $cpath = "";
my $usage = "usage : $0 -minSize minSize -cdir out_contig_directory -cpath path < locations\n";

&GetOptions("-cdir=s" => \$dirC,
            "-minSize=i" => \$minSize,
            "-cpath=s" => \$cpath,
            );

(length($cpath) > 0 && length($dirC) > 0) || die $usage;

my $cmdOpts = " -min $minSize ";
my $prg = "$ENV{CRAIG_HOME}/perl/bin/splitGenes.pl";
my %feats = %{&Utils::groupByPattern(\*STDIN, '^>\S+\s+(\d<|<)?([^\[,\],;]+)_\d+_\d+')};
unlink($dirC);

foreach my $cid (keys %feats) {
    my $ctg = $cpath."/$cid.fa";
    if(!-e $ctg) {
        $ctg = $cpath."/$cid";
    }
    if(!-e $ctg) {
        warn "$ctg cannot be found\n";
        return;
    }

    my $outCtg = "$dirC.$$";
    my $myCmdOpts = $cmdOpts." -inC $ctg -outC $outCtg ";
    Utils::execProgram($prg, $feats{$cid}, $myCmdOpts, \*STDOUT);
    `fasta2files.pl -dir $dirC < $outCtg`;
    unlink($outCtg);
}
