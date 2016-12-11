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

my $fastaFile = "";
my $contigFile = "";
my $modelName = "conf";
my $usage = "usage : $0 -model model_file -ctg contigs -F fastaFile\nReads configuration from model_file and computes histogram lengths for all introns/exons in fastaFile\n";

&GetOptions("F=s" => \$fastaFile,
	    "ctg=s" => \$contigFile,
	    "model=s" => \$modelName,
	    );

($fastaFile ne "") || die "score file must be specified. ".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

use enum qw(EXON=0 INTRON=1 INTERGENIC=4);


my $model = GeneModelProps->new();
$model->retrieveProps($modelName);
my $GF = Genefinder->new('Model' => $model);
my $org = Organism->new('Name' => $model->{'Org'}, 'Domain' => $model->{'Domain'});
print STDERR "Loading sequences.. \n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(FASTA, "<$fastaFile") || die "file $fastaFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*FASTA, 'TRAIN', 1)};
my %gcClasses = %{$model->{'GCClasses'}};

my %locsByGCClass = % {$org->getGCClasses('TRAIN', \%gcClasses)};

#$GF->trainCodingRegions($org->{'Features'}->{'TRAIN'}, \%locsByGCClass);
$GF->trainIStateRegions($org, 'TRAIN');

close(CONTIG);
close(FASTA);
print STDERR "Done!\n";
