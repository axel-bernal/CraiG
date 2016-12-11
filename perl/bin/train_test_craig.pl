#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use Getopt::Long;
use strict;

my $with = "";
my $replace = "";
my $prefix_file = "";
my $usage = "$0 train_cmd test_cmd eval_cmd [-prefix file_prefix -replace replace -with with -turns turns] \n";

my $train_cmd = shift || die $usage;
my $test_cmd = shift || die $usage;
my $eval_cmd = shift || die $usage;
my $turns = "";

&GetOptions("-prefix=s" => \$prefix_file,
	    "-replace=s" => \$replace,
	    "-with=s" => \$with,
	    "-turns=s" => \$turns);

my $params_file = "params";
my $out_file = "test";
 
if(defined $prefix_file) {
    $params_file = $prefix_file.$params_file;
    $out_file = $prefix_file.$out_file;
}

open(TRAIN, $train_cmd) || die "can't open $train_cmd. $usage";
my $train_command = <TRAIN>; chomp($train_command);
close(TRAIN);

my $rv = fork();
if($rv == 0) {
    my @turns = split(/\,/, $turns);
    foreach my $i (@turns) {
	print STDERR "looking for $params_file$i.. ";
	while(!-e "$params_file$i" || -z "$params_file$i") {
	    sleep(5);
	}
	`sync`;
	sleep(5);
	print STDERR "found!\n";

	my $nparams_file = "$params_file$i";
	if($replace ne "" && $with ne "") {
	    $nparams_file = "$nparams_file.$$";
	    `sed 's/$replace/$with/g' $params_file$i > $nparams_file`;
	}

	my $testd = `sed 's/PARAMS_FILE/$nparams_file/g' $test_cmd`;
	print STDERR "executing $testd > $out_file$i.locs 2>> $out_file.log\n";
	system($testd." > $out_file$i.locs 2>> $out_file.log");
	my $evald = `sed 's/PRED_FILE/$out_file$i/g' $eval_cmd`;
	print STDERR "executing $evald\n";
	system($evald);

	if($replace ne "" && $with ne "") {
	    unlink($nparams_file);
	}
    }
    exit();
}

system("$train_command");
waitpid($rv, 0);

