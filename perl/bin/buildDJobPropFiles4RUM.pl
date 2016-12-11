#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $prefix_propdir;
my $species = "tgondii";
my $dataset_id = "oocyst_booth";
my $sample_id = "day0";
my $rsq_inputdir = "/gpfs/fs121/h/abernal/Applications/Contrib/rum/data/me49-8.0.oocyst_boothday0";
my $prefix_craig_output = "/gpfs/fs121/h/abernal/Data/me49-8.0/me49-8.0.oocyst_boothday0";
my $fasta_file = "/gpfs/fs121/h/abernal/Data/me49-8.0/TgondiiME49_ToxoDB-8.0.gff";
my $annot_file = "/gpfs/fs121/h/abernal/Data/me49-8.0/ToxoDB-8.0_TgondiiME49_Genome.fasta";

my $usage = "usage : $0 [Options]\n\t--species\tSPECIES\n\t--prefix-propdir\tPREFIX_PROPDIR\n\t--dataset-id\tDATASET_ID\n\t--sample-id\tSAMPLE_ID\n\t--rsq-inputdir\tRSQ_INPUTDIR\n\t--prefix-craig-output\tPREFIX_CRAIG_OUTPUT\n\t--fasta-file\tFASTA_FILE\n\t--annot-file\tANNOT_FILE\n";

&GetOptions("--prefix-propdir=s" => \$prefix_propdir,
	    "--species=s" => \$species,
	    "--data-set-id=s" => \$dataset_id,
	    "--sample-id=s" => \$sample_id,
	    "--rsq-inputdir=s" => \$rsq_inputdir,
	    "--prefix-craig-output=s" => \$prefix_craig_output,
	    "--fasta-file=s" => \$fasta_file,
	    "--annot-file=s" => \$annot_file,
    );

$prefix_propdir ne "" || die "PREFIX_DIR must be given\n".$usage;

my $prop_dir =  $prefix_propdir."_input";
my $djob_dir = $prefix_propdir."_DJob";

unlink($prop_dir) if -e $prop_dir;
unlink($djob_dir) if -e $djob_dir;
`mkdir $prop_dir`;
`mkdir $djob_dir`;

my $task_props = "RSqType=rum\nRSqStranded=yes\nRSqWrapperIn=cat\nInputAnnotFormat=gff3\nInputContigFormat=fasta\ntranscriptTag=exon\nCDSTag=CDS\ngcClasses=100\ncoverageDensity=2\ngeneModelType=ngscraig\nblockLen=20\nnumPermutations=1000\nmodelUTRs=yes\ngenUTROnlyModel=yes\nnumDJobNodes=40\nDJobNodeClass=SgeNode\nminPercJunctionAlignments=75\nspecies=$species\nRSqDataSetID=$dataset_id\nRSqSampleID=$sample_id\nRSqInputDir=$rsq_inputdir\ndataDirOutput=$prefix_craig_output.preproc\nmodelDirOutput=$prefix_craig_output.preproc.model\nDJobInputBaseDir=$prefix_craig_output.preproc.test\nInputAnnotFile=$fasta_file\nInputContigFile=$annot_file\n";

my $controller_props = "masterdir=$prop_dir/master\ninputdir=$prop_dir\nnodedir=$djob_dir\nslotspernode=1\nsubtasksize=1\ntaskclass=DJob::DistribJobTasks::CRAIGTask\nnodeclass=DJob::DistribJob::SgeNode\n";

open(TASK_PROP, ">$prop_dir/task.prop") || die "Can't open $prop_dir/task.prop";
print TASK_PROP $task_props;
close(TASK_PROP);
open(CONTROLLER_PROP, ">$prop_dir/controller.prop") || die "Can't open $prop_dir/controller.prop";
print CONTROLLER_PROP $controller_props;
close(CONTROLLER_PROP);
