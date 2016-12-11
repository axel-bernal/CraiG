#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Getopt::Long;
use File::Basename;
use Genefinder;
use GenUtils;
use Organism;
$| = 1;

# -*- perl -*-
my $modelName = "";
my $baseDir = "";
my $contigs = "";
my $maxNumExons = 10;
my $binDir = "$ENV{CRAIG_HOME}/Bin/C++/";
my $verbose = 1;

# Everything in fasta files
my $usage = "usage: $0 [Options]\nOptions:\n\t--verbose\tverbose mode\n\t--basedir <basename, include complete path>\n\t--model < model name>\n\t--contigs <contig file>\n";

&GetOptions(
            "contigs=s" => \$contigs,
	    "model=s" => \$modelName,
            "basedir=s" => \$baseDir,
	    "verbose" => \$verbose,
	    );

(length($modelName) != 0 && length($baseDir) != 0 && length($contigs) != 0) || die "Must specify model and contig files\n$usage";

$ENV{PATH} .=":$binDir";
my $model = GeneModelProps->new();
my $scoreDir = $baseDir;
if(!-e $scoreDir) {
    mkdir($scoreDir);
}


$model->retrieveProps("$modelName");

my $GF = Genefinder->new('Model' => $model);

my %gcClasses = %{$model->{'GCClasses'}};

open(CONTIGS,"<$contigs") || die "Cannot open $contigs file\n";
my $org = Organism->new('Name' => $model->{'ORG'}, 'Domain' => $model->{'DOMAIN'});
my %contigs = % {$org->loadContigs(\*CONTIGS)};

close(CONTIGS);

print STDERR "Done!\nAnalyzing Contig Signals... \n";
foreach my $cKey (keys %{$org->{'Contigs'}}) {
    my %contigSignals = ();
    print STDERR "Done!\nAnalyzing contig ", $cKey, "\n";
    my $contig = $org->{'Contigs'}->{$cKey};
    $contig->getSignals(\%contigSignals, $org->{'SignalPatterns'});
    $GF->setContigSignals(\%contigSignals);
    
    print STDERR "Done!\nStoring Signals...\n";
#$GF->scoreContigSignals($org, $maxNumExons, $verbose, $scoreDir);
    $GF->storeContigSigScores("$scoreDir/$cKey", $org->{'Contigs'});
}


# scoring signals and coding regions within the organism contigs                                 
print STDERR "Done!\nGenerating GC clusters...\n";
$GF->genGCIsochores($org->{'Contigs'}, \%gcClasses, $scoreDir);
#$GF->scoreGCIsochores($org->{'Contigs'}, \%gcClasses, $scoreDir);
#$GF->genContigClusters($org->{'Contigs'}, \%gcClasses, $scoreDir);
#$GF->scoreContigClusters(\%gcClasses, $scoreDir);
#$GF->storeClusters($scoreDir, $org->{'Contigs'});

# getting a copy of the contig file
open(CONTIGS_FIXED, ">$scoreDir/contigs") || die "Can't open $scoreDir/contigs\n";

&GenUtils::displayContigs($org->{'Contigs'}, \*CONTIGS_FIXED);
close(CONTIGS_FIXED);
