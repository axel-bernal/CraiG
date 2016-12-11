#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Utils;
use Getopt::Long;
$| = 1;
my $usage = "usage: $0 -check -no_alt -cpath path < gff_file > locations file\n-check attempts to correct bad annotations, such as phantom small introns (<5), incorrect exon phase, inexistent start, stop, donor or acceptor signals in the sequence.\n- no_alt option (off by default) only outputs the transcript with more exons for each gene\n-cpath path to the contig files\n";
my $cpath = "";
my $check = 0;
my $no_alt = 0;

&GetOptions("-check" => \$check,
            "-cpath=s" => \$cpath,
            "-no_alt" => \$no_alt,
            );

(length($cpath) > 0) || die $usage;

my $cmdOpts = ($check ? "-check " : "").($no_alt ? "-no_alt " : "");

my $prg = "$ENV{CRAIG_HOME}/perl/bin/gtf2Locs.pl";
my %feats = %{&Utils::groupByPattern(\*STDIN, '^(\s*)(\S+)\s+\S+')};

foreach my $cid (keys %feats) {

    my $ctg = $cpath."/$cid.fa";
    if(!-e $ctg) {
        $ctg = $cpath."/$cid";
    }
    if(!-e $ctg) {
        warn "$ctg cannot be found\n";
        exit;
    }
    my $myCmdOpts = " -ctg $ctg $cmdOpts";
    Utils::execProgram($prg, $feats{$cid}, $myCmdOpts, \*STDOUT);
}
