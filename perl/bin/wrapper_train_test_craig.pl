#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use Getopt::Long;
use strict;

#my $qsub = "nohup ";
my $qsub = "qsub -p 0 -j y -o nohup.out -q all.q -cwd -b yes -V ";
my $contigFile;
my $ds = 0;
my $ups = 0;
my $usage = "usage : $0 train_cmd test_cmd eval_cmd outdir -test test_file -opt train_opts -conf config -replace replace -with with -turns turns\n";

my $train_file = shift || die $usage;
my $test_cmd = shift || die $usage;
my $eval_cmd = shift || die $usage;
my $outdir = shift || die $usage;
my $replace = "";
my $with = "";
my $train_opt_file = "";
my $turns = "";
my $conf = "";
my $test = "";

&GetOptions("-conf=s" => \$conf,
	    "-test=s" => \$test,
            "-opt=s" => \$train_opt_file,
	    "-replace=s" => \$replace,
	    "-with=s" => \$with,
	    "-turns=s" => \$turns,
    );


open(TRAIN, $train_file) || die "can't open $train_file. $usage";
my $line;
my $train_cmd = "";
while(defined ($line = <TRAIN>)) {
    if($line =~ /^#/) {
	next;
    }
    $train_cmd = $train_cmd.$line;
    chomp($train_cmd);
}
close(TRAIN);

my @train_opts = ();
print STDERR $train_opt_file, "\n";
if($train_opt_file =~ /\S+/ && $train_opt_file ne "none") {
    @train_opts = `cat $train_opt_file`;
}
else {
    push(@train_opts, "");
}

my $addoptions = " -turns $turns ";

if($replace ne "" && $with ne "") {
    $addoptions = " -replace $replace -with $with".$addoptions;
}

for(my $i = 0; $i < scalar(@train_opts); $i++) {
    my $train_opt = $train_opts[$i];
    chomp($train_opt);
    my $od = "$outdir";
    my $mytest = $test;

    if(scalar(@train_opts) > 1) {
	$mytest =~ s/ITERATION/$i/g;
	$od = "$outdir.$i";
    }

    my $myaddoptions = $addoptions;
    if($mytest ne "") {
	$myaddoptions = "-test $mytest".$addoptions;
    }
    
    if(!-e "$od") {
	`mkdir "$od"`;
    }
    else {
	system("rm $od/*");
    } 

    my $command = "nohup craigTrain $train_opt $train_cmd";
    $command =~ s/ITERATION/$i/g;

    if($conf ne "") {
	$command =~ s/CONF_FILE/$conf/g;
    }
    
    my $prg_file = "$train_file.$$.$i";
    open(TMP, ">$od/$prg_file");
    print TMP $command;
    close(TMP);
    `chmod +x $od/$prg_file`;
    print STDERR "cd $od && $qsub $prg_file\n";
    print STDERR "cd $od && nohup test_eval_craig.pl ../$test_cmd ../$eval_cmd -sleep 5 $myaddoptions 2> ../$od/test.log\n";
    `cd $od && $qsub $prg_file`;
    `cd $od && nohup test_eval_craig.pl ../$test_cmd ../$eval_cmd -sleep 5 $myaddoptions 2>> ../$od/test.log`;
    unlink($prg_file);
}

