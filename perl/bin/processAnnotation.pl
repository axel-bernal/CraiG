#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use GeneModelProps;

my $usage = "$0 [-s contigs locs] [-t test_set] [gc_classes]]\n pred_file must have a .locs suffix and must given with absolute path. test_set is any of {1-GP1,2-GP2,3-BG,4-TIGR,5-GEN178,6-H178,7-HMR, 8-TEST65, 9-MyGEN178, 10-COMPTEST65, 11-BG_HMR_EMBL, 12-Encode_test, 13-Encode, 14-GP1+GP2, 15-ShortIntrons, 16-test, 17-comptest}\ngc_classes is \"43,51,57\" by default";

my $tag = shift || die "$usage\n";
my $contigs;
my $fasta;
my $GCClasses = "43,51,57";
my $t = -1;
my $pwd = `pwd`;
chomp $pwd;
if($tag eq "-s") {
    $contigs = shift || die $usage;
    $fasta = shift || die $usage;
    $contigs = $pwd."/".$contigs;
    $fasta = $pwd."/".$fasta;
}    
elsif($tag eq "-t") {
    $t = shift || die;
}
 
if($#ARGV >= 0) { $GCClasses = shift; }

@args = ();
@{$args[1]} = ("~/Output/HUMAN1/geneparser_set_I", "~/Data/Human/GeneParser/test_set_I", "~/Data/Human/GeneParser/test_set_I");
@{$args[2]} = ("~/Output/HUMAN1/geneparser_set_II", "~/Data/Human/GeneParser/test_set_II", "~/Data/Human/GeneParser/test_set_II");
@{$args[3]} = ("~/Output/HUMAN1/bursguig", "~/Data/Human/BursetGuigo/BursetGuigo", "~/Data/Human/BursetGuigo/BursetGuigo");
@{$args[4]} = ("~/Output/HUMAN1/tigr", "~/Data/Human/Tigr/Test/tigr_test_cannon", "~/Data/Human/Tigr/Test/tigr_test");
@{$args[5]} = ("~/Output/HUMAN1/Gen178", "~/Data/Human/Guigo2000/Gen178", "~/Data/Human/Guigo2000/Gen178");
@{$args[6]} = ("~/Output/HUMAN1/embl50.h178", "~/Data/Human/Guigo2000/embl50.h178", "~/Data/Human/Guigo2000/embl50.h178");
@{$args[7]} = ("~/Output/HUMAN1/hmr", "~/Data/Human/HMR/hmrgen", "~/Data/Human/HMR/hmrgen");
@{$args[8]} = ("~/Output/HUMAN1/test_65", "~/Data/Human/GenScan/genscan_testing", "~/Data/Human/GenScan/genscan_testing");
@{$args[9]} = ("~/Output/HUMAN1/MyGEN178", "~/Data/Human/MyGEN178/MyGEN178", "~/Data/Human/MyGEN178/MyGEN178");
@{$args[10]} = ("~/Output/HUMAN1/comptest_65", "~/Data/Human/GenScan/genscan_testing_comp", "~/Data/Human/GenScan/genscan_testing_comp");
@{$args[11]} = ("~/Output/HUMAN1/comptest_65", "~/Data/Human/BG_HMR_EMBL/bg_hmr_embl", "~/Data/Human/BG_HMR_EMBL/bg_hmr_embl");
@{$args[12]} = ("", "~/Data/Human/Encode/programs/perceptron/encode_rstest", "~/Data/Human/Encode/encode_test");
@{$args[13]} = ("", "~/Data/Human/Chromo/Vega/Encode/encode_test", "~/Data/Human/Chromo/Vega/Encode/encode_test");
#@{$args[13]} = ("", "~/Data/Human/Encode/encode", "~/Data/Human/Encode/encode");
@{$args[14]} = ("", "~/Data/Human/GP1+GP2/test", "~/Data/Human/GP1+GP2/test");
@{$args[15]} = ("", "~/Data/Human/ShortIntronGenes/shortIntrons", "~/Data/Human/ShortIntronGenes/shortIntrons");
@{$args[16]} = ("", "~/Data/Human/Encode/encode_test_noalt_trimmed", "~/Data/Human/Encode/encode_test_noalt_trimmed");
@{$args[17]} = ("", "~/Data/Human/Encode/encode_test_alt_trimmed", "~/Data/Human/Encode/encode_test_alt_trimmed");
@{$args[18]} = ("", "~/Data/Human/Tigr/SignalPeptide/tigr_nopep", "~/Data/Human/Tigr/SignalPeptide/tigr_nopep");
@{$args[19]} = ("", "~/Data/Human/Tigr/SignalPeptide/tigr_pep", "~/Data/Human/Tigr/SignalPeptide/tigr_pep");
@{$args[20]} = ("", "~/Data/Human/All/train_genes", "~/Data/Human/All/train_contigs.xsmall");
@{$args[21]} = ("", "~/Data/Celegans/nGASP/Validation/genes.nalt", "~/Data/Celegans/nGASP/Validation/contigs.nalt.xsmall.masked");
@{$args[22]} = ("", "~/Data/Celegans/nGASP/Validation/genes.nalt.trimmed", "~/Data/Celegans/nGASP/Validation/contigs.nalt.trimmed.xsmall.masked");
@{$args[23]} = ("", "~/Data/Celegans/GAZE/WormSeq", "~/Data/Celegans/GAZE/WormSeq.xsmall.masked");
@{$args[24]} = ("", "~/Data/Celegans/GAZE/WormSeq.split", "~/Data/Celegans/GAZE/WormSeq.split.xsmall.masked");
@{$args[25]} = ("", "~/Data/Celegans/GAZE/WormSeq.split2", "~/Data/Celegans/GAZE/WormSeq.split2.xsmall.masked");
@{$args[26]} = ("", "~/Data/Celegans/nGASP/Test/genes", "~/Data/Celegans/nGASP/Test/test_regions.xsmall.masked");
@{$args[27]} = ("", "~/Data/Celegans/nGASP/Test/allgenes", "~/Data/Celegans/nGASP/Test/test_regions.xsmall.masked");
@{$args[28]} = ("", "~/Data/Human/Chromo/Vega/naltraincodgenes.trimmed", "~/Data/Human/Chromo/Vega/naltraincontigs.trimmed");
@{$args[29]} = ("", "~/Data/Human/AllIGenic/train_genes", "~/Data/Human/AllIGenic/train_contigs.xsmall");
@{$args[30]} = ("", "~/Data/Human/Chromo/Vega/naltraincodgenes.nr.trimmed", "~/Data/Human/Chromo/Vega/naltraincontigs.nr.trimmed");
@{$args[31]} = ("", "~/Data/Human/Augustus/human.train", "~/Data/Human/Augustus/human.train");
@{$args[32]} = ("", "~/Data/Human/GenScan/genscan_training", "~/Data/Human/GenScan/genscan_training.xsmall");
@{$args[33]} = ("", "~/Data/Human/Chromo/Vega/Uniprot/NSPP.run.0/NSPP.test", "~/Data/Human/Chromo/Vega/Uniprot/NSPP.run.0/NSPP.test");
@{$args[34]} = ("", "~/Data/Human/Chromo/Vega/Encode/encode_train.geo", "~/Data/Human/Chromo/Vega/Encode/encode_train.geo");
@{$args[35]} = ("", "~/Data/Celegans/Peptide/NSPP.run.0/NSPP.valid", "~/Data/Celegans/Peptide/NSPP.run.0/NSPP.valid");
@{$args[36]} = ("", "~/Data/Celegans/ConfirmedGenes/nGASP/test_codgenes.trimmed", "~/Data/Celegans/ConfirmedGenes/nGASP/test_contigs.trimmed");
@{$args[37]} = ("", "~/Data/Human/Chromo/Vega/test.run.0/codgenes.nov.test", "~/Data/Human/Chromo/Vega/test.run.0/contigs.nov.test");
@{$args[38]} = ("", "~/Data/toxo43/tgondii.test", "~/Data/toxo43/tgondii.test");
@{$args[39]} = ("", "~/Data/toxo43/tgondii.final.run.0/tgondii.final.valid", "~/Data/toxo43/tgondii.final.run.0/tgondii.final.valid");
@{$args[40]} = ("", "~/Data/toxo43/compare/tgondii.test.annot", "~/Data/toxo43/tgondii.test", );
@{$args[41]} = ("", "~/Data/toxo43/tgondii.test", "~/Data/toxo43/tgondii.test", );
@{$args[42]} = ("", "~/Data/toxo43/tgondii.test.perfect", "~/Data/toxo43/tgondii.test.perfect", );
@{$args[43]} = ("", "~/Data/plasmo55/plasmo.test.complete", "~/Data/plasmo55/plasmo.test.complete.masked", );
@{$args[44]} = ("", "~/Data/plasmo55/plasmo.tests/plasmo.test.complete.0", "~/Data/plasmo55/plasmo.tests/plasmo.test.complete.0", );
@{$args[45]} = ("", "~/Data/plasmo55/plasmo.tests/plasmo.test.complete.1", "~/Data/plasmo55/plasmo.tests/plasmo.test.complete.1", );
@{$args[46]} = ("", "~/Data/Athaliana/Combiner/test-genes", "~/Data/Athaliana/Combiner/test-genes", );
@{$args[47]} = ("", "~/Data/Athaliana/Combiner/train-genes.trim", "~/Data/Athaliana/Combiner/train-genes.trim", );
@{$args[48]} = ("", "~/Data/Athaliana/Augustus/araset", "~/Data/Athaliana/Augustus/araset", );

if($tag ne "-s") {
    $contigs = "$args[$t]->[2]";
    $fasta = "$args[$t]->[1]";
}

my $gtf = "$fasta.gtf";
my $model = GeneModelProps->new();
my %gcClasses = %{$model->setGCClasses($GCClasses)};

system("$ENV{CRAIG_HOME}/perl/bin/divideLocsByGCClass.pl -createOffset -ctg $contigs.fa -GC $GCClasses -F $fasta.locs");

if(-e "$gtf.modcoord") { 
    system("rm $gtf.modcoord");
}

my $gcClass;
while(($gcClass, undef) = each (%gcClasses)) {
    system("$ENV{CRAIG_HOME}/perl/bin/locs2GTF.pl $contigs.fa < $fasta.locs.$gcClass > $gtf.$gcClass");
    system("cat $gtf.$gcClass >> $gtf.modcoord");
    system("$ENV{CRAIG_HOME}/perl/bin/get_general_stats.pl -g -q $gtf.$gcClass > $fasta.stats.$gcClass");
}
system("$ENV{CRAIG_HOME}/perl/bin/get_general_stats.pl -g -q $gtf.modcoord > $fasta.stats");
