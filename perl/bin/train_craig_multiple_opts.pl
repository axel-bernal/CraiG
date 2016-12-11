#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use strict;

my $usage = "$0 train_cmd train_opts dir_p\n";

my $train_file = shift || die $usage;
my $train_opt_file = shift || die $usage;
my $dir_p = shift || die $usage;

open(TRAIN, $train_file) || die "can't open $train_file. $usage";
my $train_cmd = <TRAIN>; chomp($train_cmd);
close(TRAIN);
my @train_opts = `cat $train_opt_file`;

for(my $i = 0; $i < scalar(@train_opts); $i++) {
    my $train_opt = $train_opts[$i];
    chomp($train_opt);

    if(!-e "$dir_p.$i") {
	`mkdir "$dir_p.$i"`;
    }
    else {
	`rm "$dir_p.$i/*"`;
    } 
    
    my $command = "nohup craigTrain $train_opt $train_cmd";
    open(TRAIN, ">$train_file.$i") || die "$train_file.$i\n";
    print TRAIN $command;
    close(TRAIN);
    print STDERR "cd $dir_p.$i && $command > nohup.out\n";
    `cd $dir_p.$i && nohup train_test_craig.pl ../$train_file.$i ../test_cmd ../eval_cmd 2> ../$dir_p.$i/nohup.out`;
    unlink("$train_file.$i");
}
