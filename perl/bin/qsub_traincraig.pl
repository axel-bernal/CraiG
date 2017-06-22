#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use Getopt::Long;
use strict;

my $usage = "$0 conf --best=topK --mpi=<np> --parallel-comb=euclid|k-l --comb=euclid|k-l --avg=none|last|all|weighted --ml=MIN|EXP|SEP|FIRST --orw=ocw --oc=TOP|BETTER|ALL --dev --tm=trainin_method --r=r --phi=phi --limit=limit --p=priority --rs=restart_file --loss=<loss_function> --bs=<blockSize> -lb\n";

my $conf = shift || die $usage;
my $tm = "MIRA";
my $balanceLoad;
my $phi = "0.1";
my $avg = "all";
my $comb_multi = "euclid";
my $comb = "euclid";
my $r = "1.0";
my $orw = "-1";
my $oc = "NONE";
my $limit = "15";
my $priority = "0";
my $mpi = 1;
my $restart = "";
my $dev = 0;
my $best = "1";
my $addfeats;
my $multilab = "ALL";
my $loss = "HAMMING";
my $blockSize = "1";
&GetOptions("-tm=s" => \$tm,
            "-dev" => \$dev,
            "-a" => \$addfeats,
            "-orw=s" => \$orw,
            "-oc=s" => \$oc,
            "-bl" => \$balanceLoad,
            "-bs=s" => \$blockSize,
            "-best=s" => \$best,
            "-phi=s" => \$phi,
            "-mpi=i" => \$mpi,
            "-r=s" => \$r,
            "-avg=s" => \$avg,
	    "-comb=s" => \$comb,
	    "-parallel-comb=s" => \$parallel_comb,
            "-limit=s" => \$limit,
            "-p=s" => \$priority,
            "-rs=s" => \$restart,
            "-ml=s" => \$multilab,
            "-loss=s" => \$loss,
            );

my $output = "nohup.out";
if(length($restart) > 0) {
    $restart = " --restart=$restart ";
}


if($addfeats) {
    $addfeats = "-a";
}
else {
    $addfeats = "";
}

if($balanceLoad) {
    $balanceLoad = "-bl";
}
else {
    $balanceLoad = "";
}

my $stringMPI = "";
if($mpi > 1) {
    $stringMPI = "-pe OpenMPI $mpi mpirun -np $mpi";
}

my $command;

if($dev) {
    $command = "qsub -p $priority -j y -o $output -q all.q -cwd -b yes -V $stringMPI $ENV{HOME}/Projects/newmpicraig/bin/craigTrain --block-size=$blockSize --oracle-upd=$oc --oracle-wait=$orw --best=$best --multi-label-upd=$multilab --r=$r -v $balanceLoad $addfeats $restart --algorithm=$tm --phi=$phi --limit=$limit --avg-method=$avg --parallel-comb-method=$parallel_comb --comb-method=$comb --loss-function=$loss $conf params";
} 
else {
    $command = "qsub -p $priority -j y -o $output -q all.q -cwd -b yes -V $stringMPI $ENV{CRAIG_HOME}/bin/craigTrain --block-size=$blockSize --oracle-upd=$oc --oracle-wait=$orw --best=$best --multi-label-upd=$multilab --r=$r -v $balanceLoad $addfeats $restart --algorithm=$tm --phi=$phi --limit=$limit --avg-methodv=$avg --parallel-comb-method=$parallel_comb --comb-method=$comb --loss-function=$loss $conf params";
}
print $command, "\n";
system($command);
