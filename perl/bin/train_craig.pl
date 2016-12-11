#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use Getopt::Long;
use strict;

my $qsub = "nohup ";
my $qsub = "qsub -p 0 -j y -o nohup.out -q all.q -cwd -b yes -V ";
my $contigFile;
my $ds = 0;
my $ups = 0;
my $usage = "usage : $0 train_cmd outdir -opt train_opts -conf config\n";

my $train_file = shift || die $usage;
my $outdir = shift || die $usage;
my $train_opt_file = "";
my $conf = "";

&GetOptions("-conf=s" => \$conf,
            "-dir=s" => \$outdir,
            "-opt=s" => \$train_opt_file,
    );


open(TRAIN, $train_file) || die "can't open $train_file. $usage";
my $train_cmd = <TRAIN>; chomp($train_cmd);
close(TRAIN);

my @train_opts = ();
print STDERR $train_opt_file, "\n";
if($train_opt_file =~ /\S+/ && $train_opt_file ne "none") {
    @train_opts = `cat $train_opt_file`;
}
else {
    push(@train_opts, "");
}


for(my $i = 0; $i < scalar(@train_opts); $i++) {
    my $train_opt = $train_opts[$i];
    chomp($train_opt);
    my $od = "$outdir";

    if(scalar(@train_opts) > 1) {
	$od = "$outdir.$i";
    }

    if(!-e "$od") {
	`mkdir "$od"`;
    }
    else {
	system("rm $od/*");
    } 

    my $command = "craigTrain $train_opt $train_cmd";
    if($conf ne "") {
	$command =~ s/CONF_FILE/$conf/g;
    }
    
    my $prg_file = "$train_file.$$.$i";
    open(TMP, ">$od/$prg_file");
    print TMP $command;
    close(TMP);
    `chmod +x $od/$prg_file`;
    print STDERR "proc# $$: cd $od && $command > nohup.out\n";
    `cd $od && $qsub $prg_file`;    
#    unlink($prg_file);
}
