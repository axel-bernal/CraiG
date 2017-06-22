#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use strict;

my $usage = "$0 train_cmd train_opts\n";

my $train_file = shift || die $usage;
my $train_opt_file = shift || die $usage;

open(TRAIN, $train_file) || die "can't open $train_file. $usage";
my $train_cmd = <TRAIN>; chomp($train_cmd);
close(TRAIN);
my @train_opts = `cat $train_opt_file`;

for(my $i = 0; $i < scalar(@train_opts); $i++) {
    my $train_opt = $train_opts[$i];
    chomp($train_opt);

    if(!-e "arrow.$offset") {
	`mkdir "arrow.$offset"`;
    }
    else {
	`rm "arrow.$offset/*"`;
    } 
    
    my $command = "qsub -p 0 -j y -o nohup.out -q all.q -cwd -b yes -V $stringMPI $ENV{CRAIG_HOME}/bin/train_craig $train_opt $train_cmd";

    print STDERR "cd arrow.$offset && $command\n";
    system("cd arrow.$offset && $new_train_cmd");
}
