#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use Getopt::Long;
use strict;

my $out_file = "test";
#my $qsub = "nohup ";
my $qsub = "qsub -p 0 -j y -o $out_file.log -q all.q -cwd -b yes -V ";
my $with = "";
my $replace = "";
my $prefix_file = "";
my $sleep = 0;
my $usage = "$0 test_cmd eval_cmd [-sleep sleep_time -test test_file -prefix file_prefix -replace replace -with with -turns turns] \n";

my $test_cmd = shift || die $usage;
my $eval_cmd = shift || die $usage;
my $turns = "";
my $test = "";

&GetOptions("-prefix=s" => \$prefix_file,
	    "-sleep=i" => \$sleep,
	    "-test=s" => \$test,
	    "-replace=s" => \$replace,
	    "-with=s" => \$with,
	    "-turns=s" => \$turns);

my $test_string;
my $eval_string;

if($test ne "") {
    $eval_string = `cat $eval_cmd`;
    chomp $eval_string;
    $eval_string =~ s/EVAL_FILE/$test/g; 
    $test_string = `cat $test_cmd`;
    chomp $test_string;
    $test_string =~ s/FASTA_FILE/$test/g;
}

my @turns = ();
if($turns =~ /all(\d+)/) {
    for(my $j = 1; $j <= $1; $j++) {
	push(@turns, $j);
    }    
}
elsif($turns =~ /all/) {
    for(my $j = 1; $j <= 10; $j++) {
	push(@turns, $j);
    }
}
elsif($turns =~ /(\d+)\-(\d+)/) {
    for(my $j = $1; $j <= $2; $j++) {
	push(@turns, $j);
    }
}
else {
    @turns = split(/\,/, $turns);
}

my $params_file = "params";
 
if(defined $prefix_file) {
    $params_file = $prefix_file.$params_file;
    $out_file = $prefix_file.$out_file;
}

foreach my $i (@turns) {
    my $testd = $test_string;
    my $evald = $eval_string;

    if(!$sleep) {
	print STDERR "looking for $params_file$i.. ";
    }
    while(!-e "$params_file$i" || -z "$params_file$i") {
	sleep($sleep);
    }
    `sync`;
    sleep($sleep);
    print STDERR "found!\n";
    
    my $nparams_file = "$params_file$i";
    if($replace ne "" && $with ne "") {
	$nparams_file = "$nparams_file.$$";
	print STDERR "sed $replace $with";
	`sed 's/$replace/$with/g' $params_file$i > $nparams_file`;
    }
    
    $testd =~ s/PARAMS_FILE/$nparams_file/g;
    $testd =~ s/PRED_FILE/$out_file$i/g;
    print STDERR "executing $testd 2>> $out_file.log\n";
    $evald =~ s/PRED_FILE/$out_file$i/g;
    print STDERR "executing $evald\n";
    
    my $prg_file = "test_cmd.$$.$i";
    open(TMP, ">$prg_file");
    print TMP $testd,";\n";
    print TMP $evald;
    close(TMP);
    `chmod +x $prg_file`;
    `$qsub $prg_file`;
#    unlink($prg_file);
    
    if($replace ne "" && $with ne "") {
#	unlink($nparams_file);
    }
}


