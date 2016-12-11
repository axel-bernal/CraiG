#!/usr/bin/perl
###############################################################################
# train_model.pl by Axel Bernal
###############################################################################

($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Getopt::Long;
use File::Basename;
use Utils;
use Genefinder;
use GeneModelProps;
use Organism;
$| = 1;

use enum qw(START=0 STOP=1 DONOR=2 ACEPTOR=3);
use enum qw(STRAND_FWD STRAND_COMP);

# -*- perl -*-
my $model_name = "";
my $fmtType = "fasta";
my $genes = "";
my $contigs = "";
my $evalGenes = "";
my $evalContigs = "";
my $starts = undef;
my $donors = undef;
my $aceptors = undef;
my $stops = undef;
my ($startWindow, $stopWindow, $donorWindow, $aceptorWindow) = ("7:8", "3:18","3:6", "4:5");
# defaults for parameters
my $orgName = "HS";
my $domain = "E";
my $GCClasses = "43,51,57";
my $binDir = "$ENV{CRAIG_HOME}/Bin/C++/";
my $includeSeq = 0;
my $maxNumExons = 100;
my $verbose = 1;
my $addCoding = "";
my $margin = -1;
my $class;

# Everything in fasta files
my $usage = "usage: $0 [Options]\nOptions:\n\t--org <orgName>\n\t--domain\t<domain>\n\t--genes\t<genes>\n\t--includeSeq\tSpecify if gene sequences are present in genes file\n\t--contigs\t<contig_file>\n\t--model <directory>\tdirectory where to find the model files\n\t--GC\t<GC delimiters>. ex:'41,50,60'\n\t--wsrt\t<window_start> ex '20-0'\n\t--wacp\t<window_aceptor>\n\t--wdon\t<window_donor>\n\t--wstp\t<window_stop>\n\t--binDir <directory> Specifies the directory containing the executable programs to be used by the genefinder\n\t--addCoding\t<Additional coding sequences>\n\t--evalContigs\t<crossEvaluation Contigs>\n\t--evalGenes\t<crossEvaluation Gene Set>\n\t--margin\t<margin>\n\t--fmt\t<gtf|fasta|genbank>\n\t--verbose\tverbose mode\n";

&GetOptions(
            "org=s" => \$orgName,
	    "domain=s" => \$domain,
            "genes=s" => \$genes,
	    "contigs=s" => \$contigs,
            "evalGenes=s" => \$evalGenes,
            "evalContigs=s" => \$evalContigs,
            "model=s" => \$model_name,
	    "GC=s" => \$GCClasses,
	    "wsrt=s" => \$startWindow,
	    "wacp=s" => \$aceptorWindow,
	    "wdon=s" => \$donorWindow,
	    "wstp=s" => \$stopWindow,
	    "includeSeq" => \$includeSeq,
	    "verbose" => \$verbose,
	    "addCoding=s" => \$addCoding,
	    "fmt=s" => \$fmtType,
	    "margin=i" =>\$margin,
	    );

(length($genes) != 0 && length($model_name) != 0 && length($contigs) != 0 && length($evalGenes) != 0 && length($evalContigs) != 0) || die "Must specify gene list and model files\n".$usage;

my $scoreFile = $model_name."/TrainContigScores";
my $model = GeneModelProps->new(
				'Name' => $model_name,
				'Org' => $orgName,
				'Domain' => $domain,
				'PromoterWindow' => "200:0",
				'SignalWindows' => [$startWindow, 
						    $stopWindow, 
						    $donorWindow, 
						    $aceptorWindow,
						    ],
				'SignalModels' => [$model_name."/starts.mod",
						   $model_name."/stops.mod",
						   $model_name."/donors.mod",
						   $model_name."/aceptors.mod",
						   ],
				'SignalFiles' => [$model_name."/starts",
						  $model_name."/stops",
						  $model_name."/donors",
						  $model_name."/aceptors",
						  ],
                                'IGenicGeoMeans' => [150000, 60000, 9000, 3500],
				'HistExonNum' => $model_name."/hist_exnum",
				'HistExonLen' => $model_name."/hist_exlen",
				'HistIntronLen' => $model_name."/hist_inlen",
				'HistIntergenicLen' => $model_name."/hist_interlen",
				'CodingModel' => $model_name."/immCoding",
				'IntergenicModel' => $model_name."/immIntergenic",
				'IntronModel' => $model_name."/immIntron",
				'GCProp' => $model_name."/gcp",
				'GHMMParams' => $model_name."/hmm_params",
				'SigScoreBinLen' => 0.1,
				'CodScoreBinLen' => 0.05,
				'IStScoreBinLen' => 0.05,
				'HistIntergenicBinLen' => 320,
				'HistExonBinLen' => 90,
				'HistIntronBinLen' => 240,
				'TrainScoreFile' => $scoreFile,
				'PredScoreFile' => "NULL",
				'Genes' => $model_name."/genes.fsa",
				'TrainContigs' => $model_name."/train_contigs.fsa",
				'PredContigs' => "NULL",
				'CrossEvalGenes' => $model_name."/eval_genes.fsa",
				'Margin' => $margin,
				);

my %gcClasses = %{$model->setGCClasses($GCClasses)};
$model->syncInfo();
`cat $contigs $evalContigs > "$model->{'TrainContigs'}.pre"`;
$ENV{PATH} .=":$binDir";
open(CONTIGS,"<$model->{'TrainContigs'}.pre") || die "Cannot open $contigs file\n".$usage;
open(LOCAS,"<$genes") || die "Cannot open $genes file\n".$usage;
open(XLOCAS,"<$evalGenes") || die "Cannot open $evalGenes file\n".$usage;

### Synching model. Directory structures, temporary files and so on
print STDERR "Loading Contigs and Genes... \n";

### Creating the genefinder
my $GF = Genefinder->new('Model' => $model);

my $org = Organism->new('Name' => $model->{'Org'}, 'Domain' => $model->{'Domain'});

# loading contigs and genes for training/evaluation
my $contigs = $org->loadContigs(\*CONTIGS);
my $xgenes = $org->loadGenes(\*XLOCAS, 'EVAL', $includeSeq, $fmtType);
my $genes = $org->loadGenes(\*LOCAS, 'TRAIN', $includeSeq, $fmtType);

#trim Stop Codon at the end
&GenUtils::trimStopCodon($genes);
&GenUtils::trimStopCodon($xgenes);


#storing genes and contigs with fixed-to-zero coordinates
print STDERR "Done!\nStoring properties/files... \n";
open(GENES_TRAIN, ">$model->{'Genes'}") || die "Can't open $model->{'Genes'}\n";
open(GENES_EVAL, ">$model->{'CrossEvalGenes'}") || die "Can't open $model->{'CrossEvalGenes'}\n";
open(CONTIGS_FIXED, ">$model->{'TrainContigs'}") || die "Can't open $model->{'TrainContigs'}\n";

&GenUtils::displayGeneList($org->{'Features'}->{'EVAL'}, \*GENES_EVAL, $org->{'Contigs'}, 1);
&GenUtils::displayGeneList($org->{'Features'}->{'TRAIN'}, \*GENES_TRAIN, $org->{'Contigs'}, 1);
&GenUtils::displayContigs($org->{'Contigs'}, \*CONTIGS_FIXED);
print STDERR "Done!\n";
`rm "$model->{'TrainContigs'}.pre"`;

#storing properties
$model->storeProps("$model->{'Name'}/conf");
close(GENES_TRAIN);
close(GENES_EVAL);
close(CONTIGS_FIXED);

my %locsByGCClass = % {$org->getGCClasses('TRAIN', \%gcClasses)};

### closing files
close(CONTIGS);
close(LOCAS);
close(XLOCAS);

print STDERR "Building Histograms...\n";

#building histograms
my($hstExonLen, $hstNumExons, $hstIntronLen, $hstIntergenicLen) = @ {$GF->computeHistograms($org, 'TRAIN', \%locsByGCClass)};

foreach $class (keys %locsByGCClass) {
    &Utils::hashToFile($hstIntronLen->{$class}->{'SmoothValues'}, $model->{'HistIntronLen'}.".$class", " ");
    &Utils::hashToFile($hstIntergenicLen->{$class}->{'SmoothValues'}, $model->{'HistIntergenicLen'}.".$class", " ");
    &Utils::hashToFile($hstNumExons->{'SmoothValues'}, $model->{'HistExonNum'}.".$class"," ");
    &Utils::hashToFile($hstExonLen->{'SmoothValues'}, $model->{'HistExonLen'}.".$class", " ");
}

`rm "$scoreFile.signals"`;
print STDERR "Done!\nAnalyzing Contig Signals... \n";
foreach my $cKey (keys %{$org->{'Contigs'}}) {
    my %contigSignals = ();
    print STDERR "Done!\nAnalyzing contig ", $cKey, "\n";
    my $contig = $org->{'Contigs'}->{$cKey};
    $contig->getSignals(\%contigSignals, $org->{'SignalPatterns'});
    $GF->setContigSignals(\%contigSignals);
                                                           
    print STDERR "Done!\nProcessing Signal Information...\n";
    $GF->classifySignalsWGeneSet($org, 'TRAIN');
    $GF->classifySignalsWGeneSet($org, 'EVAL');
    print STDERR "Done!\n";                                                                                
#training the signals. Output will be put in model/signal
#    print STDERR "Training Signals...\n";
#$GF->trainContigSignals($org, \@addSignals, 40, STRAND_FWD);
#    print STDERR "Done!\nStoring Signals...\n";
#    $GF->scoreContigSignals($org, $maxNumExons, $verbose, $scoreFile);
    $GF->storeContigSigScores("$scoreFile.$cKey", $org->{'Contigs'});
    `cat "$scoreFile.$cKey.signals" >> "$scoreFile.signals"`;
    `rm "$scoreFile.$cKey.signals"`;
}


#training the coding_regions
#print STDERR "Training Coding Regions...\n";
#if($addCoding ne "") {
#    open(ADD_CODING, "<$addCoding") || die "$addCoding\n";
#    my $codGenes = $org->loadGenes(\*ADD_CODING, 'TRAIN', 1);
#    &GenUtils::trimStopCodon($codGenes);
#    close(ADD_CODING);
#    %locsByGCClass = % {$org->getGCClasses('TRAIN', \%gcClasses)};
#}

#$GF->trainCodingRegions($org->{'Features'}->{'TRAIN'}, \%locsByGCClass);
#$GF->trainIStateRegions($org, 'TRAIN');
#print STDERR "Done!\n";

# scoring signals and coding regions within the organism contigs
# evaluation set

print STDERR "Done!\nScoring contigs regions...\n";
# training set
#$GF->genContigClusters($org->{'Contigs'}, \%gcClasses, $scoreFile);
#$GF->storeClusters($scoreFile, $org->{'Contigs'});

#$GF->scoreContigClusters(\%gcClasses, $scoreFile);

print STDERR "Done!\n";
#Calling training module for each subset
