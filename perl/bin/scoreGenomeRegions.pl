#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Gene;
use GeneModelProps;
use Organism;
use Contig;
use Genefinder;
use Getopt::Long;
use strict;

my $contigFile = "";
my $modelName = "conf";
my $usage = "usage : $0 -model model_file -ctg contigs\nReads configuration from model_file and scores all regions in contigs class by class\n";

&GetOptions("ctg=s" => \$contigFile,
	    "model=s" => \$modelName,
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

use enum qw(EXON=0 INTRON=1 INTERGENIC=4);


my $model = GeneModelProps->new();
$model->retrieveProps($modelName);
my $GF = Genefinder->new('Model' => $model);
my $org = Organism->new('Name' => $model->{'Org'}, 'Domain' => $model->{'Domain'});
print STDERR "Loading sequences.. \n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %gcClasses = %{$model->{'GCClasses'}};

$GF->genContigClusters($org->{'Contigs'}, \%gcClasses, $model->{'TrainScoreFile'});
$GF->storeClusters($model->{'TrainScoreFile'}, $org->{'Contigs'}); 
$GF->scoreContigClusters(\%gcClasses, $model->{'TrainScoreFile'});

close(CONTIG);
print STDERR "Done!\n";
