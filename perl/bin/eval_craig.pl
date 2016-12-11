#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use Getopt::Long;
use strict;

my $eval_cmd = "";
my $usage = "$0 annotation[no ext] prediction[no ext] -cmd eval_cmd\n";
my $test_file = shift || die $usage;
my $pred_file = shift || die $usage;
&GetOptions("-cmd=s" => \$eval_cmd);

my $evald = "processPrediction.pl PRED_FILE -s EVAL_FILE EVAL_FILE 100";
if($eval_cmd ne "") {
    $evald = `cat $eval_cmd`;
    chomp $evald;
}

$evald =~ s/EVAL_FILE/$test_file/g; 
$evald =~ s/PRED_FILE/$pred_file/g;
print STDERR "executing $evald\n";
`$evald`;


