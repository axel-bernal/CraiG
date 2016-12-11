#!/usr/bin/perl -w
use Getopt::Long;
use strict;
$| = 1;
my $usage = "$0 [-t] -prg program -D conf_file -S strand -P params -f pred_file -set test_set\n pred_file must have a .fsa suffix and must given with absolute path. test_set is any of {1-GP1,2-GP2,3-BG,4-TIGR,5-GEN178,6-H178,7-HMR, 8-TEST65, 9-MyGEN178, 10-COMPTEST65}\n";  
my $program;
my $conf;
my $strand;
my $params;
my $arg;
my $t;
my $truncated = 0;
&GetOptions("-prg=s" => \$program,
	    "-t" => \$truncated,
            "-D=s" => \$conf,
	    "-f=s" => \$arg,
	    "-set=s" => \$t,
	    "-S=s" => \$strand,
	    "-P=s" => \$params,
            );

if(!defined($program) || !defined($conf)) {
    print $usage;
    exit;
}

my $pwdir = `pwd`;
chomp($pwdir);
$program = $pwdir."/".$program;
$conf = $pwdir."/".$conf;
$params = $pwdir."/".$params;
$arg = $pwdir."/".$arg;

my @args = ();
@{$args[1]} = ("~/Projects/CRAIG/Output/HUMAN1/geneparser_set_I", "~/Projects/CRAIG/Data/Human/GeneParser/test_set_I", "~/Projects/CRAIG/Data/Human/GeneParser/test_set_I");
@{$args[2]} = ("~/Projects/CRAIG/Output/HUMAN1/geneparser_set_II", "~/Projects/CRAIG/Data/Human/GeneParser/test_set_II", "~/Projects/CRAIG/Data/Human/GeneParser/test_set_II");
@{$args[3]} = ("~/Projects/CRAIG/Output/HUMAN1/bursguig", "~/Projects/CRAIG/Data/Human/BursetGuigo/BursetGuigo", "~/Projects/CRAIG/Data/Human/BursetGuigo/BursetGuigo");
@{$args[4]} = ("~/Projects/CRAIG/Output/HUMAN1/tigr", "~/Projects/CRAIG/Data/Human/Tigr/Test/tigr_test_cannon", "~/Projects/CRAIG/Data/Human/Tigr/Test/tigr_test");
@{$args[5]} = ("~/Projects/CRAIG/Output/HUMAN1/Gen178", "~/Projects/CRAIG/Data/Human/Guigo2000/Gen178", "~/Projects/CRAIG/Data/Human/Guigo2000/Gen178");
@{$args[6]} = ("~/Projects/CRAIG/Output/HUMAN1/embl50.h178", "~/Projects/CRAIG/Data/Human/Guigo2000/embl50.h178", "~/Projects/CRAIG/Data/Human/Guigo2000/embl50.h178");
@{$args[7]} = ("~/Projects/CRAIG/Output/HUMAN1/hmr", "~/Projects/CRAIG/Data/Human/HMR/hmrgen", "~/Projects/CRAIG/Data/Human/HMR/hmrgen");
@{$args[8]} = ("~/Projects/CRAIG/Output/HUMAN1/test_65", "~/Projects/CRAIG/Data/Human/GenScan/genscan_testing", "~/Projects/CRAIG/Data/Human/GenScan/genscan_testing");
@{$args[9]} = ("~/Projects/CRAIG/Output/HUMAN1/MyGEN178", "~/Projects/CRAIG/Data/Human/MyGEN178/MyGEN178", "~/Projects/CRAIG/Data/Human/MyGEN178/MyGEN178");
@{$args[10]} = ("~/Projects/CRAIG/Output/HUMAN1/comptest_65", "~/Projects/CRAIG/Data/Human/GenScan/genscan_testing_comp", "~/Projects/CRAIG/Data/Human/GenScan/genscan_testing_comp");
my $command = "$program -S $strand -D $conf -P $params -F $args[$t]->[0]";
if($truncated) {
    $command .= " -t";
}
system("$command > $arg.fsa");
system("~/Projects/CRAIG/PerlBin/divideLocsByGCClass.pl -ctg $args[$t]->[2].contigs -GC 43,51,57 -F $arg.fsa");
system("for i in `cat ~/Projects/CRAIG/Source/gc`; do ~/Projects/CRAIG/PerlBin/locs2GTF.pl $args[$t]->[2].contigs < $arg.fsa.\$i > $arg.gtf.\$i; done");
system("rm $arg.gtf");
system("for i in `cat ~/Projects/CRAIG/Source/gc`; do cat $arg.gtf.\$i >> $arg.gtf; done");
system("cd ~/Applications/Contrib/eval; for i in `cat ~/Projects/CRAIG/Source/gc`; do evaluate_gtf.pl -g -q $args[$t]->[1].gtf.\$i  $arg.gtf.\$i > $arg.eval.\$i; done");
system("cd ~/Applications/Contrib/eval; evaluate_gtf.pl -g -q $args[$t]->[1].gtf.modcoord  $arg.gtf > $arg.eval");
