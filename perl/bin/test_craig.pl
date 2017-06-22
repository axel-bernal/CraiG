#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use Getopt::Long;
use strict;

my $test_cmd = "";
my $replace = "";
my $with = "";
my $usage = "$0 params fasta[no ext] predictions[no ext] -cmd test_cmd -replace to_be_replaced -with replacement\n";
my $params_file = shift || die $usage;
my $fasta_file = shift || die $usage;
my $pred_file = shift || die $usage;
&GetOptions("-cmd=s" => \$test_cmd,
	    "-replace=s" => \$replace,
	    "-with=s" => \$with
    );

my $testd = "craigPredict --format=locs --best=1 -m PARAMS_FILE FASTA_FILE.fa > PRED_FILE.locs";
if($test_cmd ne "") {
    $testd = `cat $test_cmd`;
    chomp $testd;
}

my $nparams_file = $params_file;
if($replace ne "" && $with ne "") {
    $nparams_file = "$params_file.$$";
    print STDERR "sed $replace $with";
    `sed 's/$replace/$with/g' $params_file > $nparams_file`;
}

$testd =~ s/PARAMS_FILE/$nparams_file/g;
$testd =~ s/FASTA_FILE/$fasta_file/g;
$testd =~ s/PRED_FILE/$pred_file/g;

print STDERR "executing $testd\n";
`$testd`;


