#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use strict;

my $usage = "$0 test_cmd eval_cmd prefix_dir number [file_prefix]\n";

my $test_cmd = shift || die $usage;
my $eval_cmd = shift || die $usage;
my $dir = shift || die $usage;
my $number = shift || die $usage;
my $prefix_file = shift;
my $params_file = "params";
my $out_file = "test";
 
if(defined $prefix_file) {
    $params_file = $prefix_file.$params_file;
    $out_file = $prefix_file.$out_file;
}

for(my $i = 0; $i < $number; $i++) {
    my $d = "$dir.$i";
    my $rv = fork();
    if($rv == 0) {
	for(my $j = 2; $j < 4; $j++) {
	    print STDERR "looking for $d/$params_file$j.. ";
	    while(!-e "$d/$params_file$j" || -z "$d/$params_file$j") {
		;
	    }
	    print STDERR "found!\n";
	    sleep(5);
	    my $nparams_file = "mod.$params_file$j";
	    `sed 's/test.final/train/' $d/$params_file$j > $d/$nparams_file`;
	    my $tcmd = `sed 's/PARAMS_FILE/$nparams_file/' $test_cmd`;
	    print STDERR "executing $tcmd > $out_file$j.locs 2>> $out_file.log\n";
	    system("cd $d && $tcmd > $out_file$j.locs 2>> $out_file.log");
	    my $ecmd = `sed 's/PRED_FILE/$out_file$j/' $eval_cmd`;
	    print STDERR "executing $ecmd\n";
	    system("cd $d && $eval_cmd");
	}
	exit();
    }
    if($i >= 6) {
	waitpid($rv, 0);
    }
}
