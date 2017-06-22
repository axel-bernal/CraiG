#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use Getopt::Long;
use strict;
my $rsq_list = "";
my $base_djob_dir = "/gpfs/fs121/h/abernal/GUS/project_home/DJob/";
my $djob_dir = "";
my $base_craig_output_dir = "/gpfs/fs121/h/abernal/Data";
my $craig_output_dir = "";
my $version = "9.0";
my $fasta_file = "";
my $gff_file = "";
my $species = "tgondii";
my $usage = "$0 --species SPECIES --rsq-list RNASeq_List --djob-dir DJob_Dir(in project_home/DJob) --craig-dir Craig_OutputDir(in ~/Data) --version Version --fasta Fasta_File --gff Gff_File\n";

&GetOptions("--rsq-list=s" => \$rsq_list,
            "--djob-dir=s" => \$djob_dir,
	    "--species=s" => \$species,
            "--fasta=s" => \$fasta_file,
	    "--gff=s" => \$gff_file,
	    "--craig-dir=s" => \$craig_output_dir,
	    "--version=s" => \$version,
    );

($fasta_file ne "" && $gff_file ne "" && $rsq_list ne "" && $djob_dir ne "" && $craig_output_dir ne "") || die $usage;
$djob_dir = $base_djob_dir."/".$djob_dir;
$craig_output_dir = $base_craig_output_dir."/".$craig_output_dir;

open(LIST, "$rsq_list") || die "Can't open $rsq_list\n";
open(SCRIPT, ">$djob_dir/proc_$rsq_list.sh") || die "Can't open $djob_dir/proc_$rsq_list.sh";
my $line;
my $counter = 0;
while(defined($line = <LIST>)) {
    $line =~ /(\S+)\s+(\S+)/;
    my $dbname = $1;
    my $sample = $2;
    $dbname =~ /([^\_]+)\_(\S+)/;
    my $org = $1;
    my $dataset = $2;
    my $orientation = 'N';
    $orientation = 'F' if $dataset =~ /Boothroyd/ || $dataset =~ /Sibley/ || $dataset =~ /Gregory/ || $dataset =~ /Hehl/;
    $orientation = 'NN' if $dataset =~ /Reid/;

    die "File $djob_dir does not exist\n" unless -e "$djob_dir";
    die "File $craig_output_dir/$fasta_file does not exist\n" unless -e "$craig_output_dir/$fasta_file";
    die "File $craig_output_dir/$gff_file does not exist\n" unless -e "$craig_output_dir/$gff_file";

    system("buildDJobPropFiles4Craig.pl --species $species --rsq-type rum --rsq-inputdir /gpfs/fs121/h/abernal/Applications/Contrib/rum/data/$org-$version.$dataset.$sample --rsq-orientation $orientation --dataset-id $dataset --sample-id $sample --prefix-propdir $djob_dir/$org-$version.$dataset.$sample --num-djobnodes 40 --prefix-craig-output $craig_output_dir/$org-$version.$dataset.$sample --fasta-file $craig_output_dir/$fasta_file --annot-file $craig_output_dir/$gff_file");
    my $dir = "$org-$version.$dataset.$sample"."_input";
    $dir = "../".$dir if $counter;
    print SCRIPT "cd $dir && rm -rf master && nohup distribjob --propFile controller.prop --numNodes 1 --memoryPerNode 20.0 > nohup.out;\n";
    $counter++;
}
close SCRIPT;
