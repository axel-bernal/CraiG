package DJob::DistribJobTasks::CRAIGTask;
#use lib "/home/abernal/GUS/gus_home/lib/perl";

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["RSqType", "rum", "RNASeq type rum|gsnap|epdb. Only rum currently supported"],
 ["RSqDataSetID", "", "RNASeq dataset id"],
 ["RSqSampleID", "", "RNASeq sample id"],
 ["RSqInputDir", "", "Full path to RNASeq's alignments directory"],
 ["RSqOrientation", "", "Orientation of RNASeq reads: N|F|R. Use 'N' for non-stranded RNA reads. For paired reads you need to provide with orientation for each individual read, i.e. RR, FF, FR (illumina)"],
 ["RSqWrapperIn", "cat", "Script used to read RNASeq coverage information"],
 ["InputAnnotFile", "", "Full path to genome isoform annotation"],
 ["InputAnnotFormat","gff3","Format for genome isoform annotation locs|gff3|gtf|genbank"],
 ["InputContigFile","","Full path to genome contig sequences"],
 ["InputContigFormat","fasta","Format for genome isoform annotation fasta|genbank"],
 ["closestSpecies","none", "Closest-related organism that has been defined in the models directory. Needs to be given if the current species has no entry in said directory"],
 ["transcriptTag", "exon", "Tag for recognizing transcribed regions in InputAnnotFile"],
 ["CDSTag", "CDS", "Tag for recognizing coding regions in InputAnnotFile"],
 ["gcClasses", "100", "gc isochores present in the genome. f.e. human would be 43,51,57"],
 ["species", "", "name of species tgondii|plasmo"],
# ["coverageDensity", "-1", "Sample's coverage density. A number between 1(very sparse, like gregory datasets) to 9(very dense, like reid day4 dataset). This is roughyl a linear factor, but most datasets with good mapping should always be above 4 or 5. For example, reid3 should around 7, sibley around 5, oocyst(the bottom line for good datasets) around 4. Gregory data sets should all be 1 or 2. This can be extracted from the mapping_stats.txt file in rum. Just check the number of total mapped reads, consider paired reads to be better than single reads and land in a number. This estimated number is not so much important as long as it's a relatively good guess. If absent the program will try to estimate one by default, using the total number of junctions that have been covered"],
 ["geneModelType", "ngscraig", "model type ngscraig/ecraig/craig"],
 ["trainingIterations", "10", "number of training iterations for learning convergence"],
 ["modelUTRs", "yes", "models and predicts UTRs as part of genes"],
 ["prefixOutputDir", "", "prefix for the directories with learnt and predicted models"],
 ["genUTROnlyModel", "yes", "generates utr_only predictions"],
 ["numDJobNodes", "20", "Number of nodes used by internal distribjob calls"],
 ["DJobNodeClass", "LocalNode", "Node Class for the internal DJob calls, either LocalNode or SgeNode"],
 ["minPercJunctionAlignments", "", "Minimum percentage of genes with introns that have good junction alignments that are needed in order to use the junction information as the sole source for splicing"]
);

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;
  # creating preconfig file
    open(PRECONFIG,">$inputDir/preconfig");
    my $rsqType = $self->getProperty('RSqType');
    my $rsqDataSetId = $self->getProperty('RSqDataSetID');
    my $rsqSampleId = $self->getProperty('RSqSampleID');
    my $rsqInputDir = $self->getProperty('RSqInputDir');
    my $rsqOrientation = $self->getProperty('RSqOrientation');
    my $rsqWrapperIn = $self->getProperty('RSqWrapperIn');
    print PRECONFIG "rnaseq\t$rsqType\t$rsqDataSetId\t$rsqSampleId\t$rsqInputDir\t$rsqOrientation\t$rsqWrapperIn\n";
    close(PRECONFIG);
    
    my $inputAnnotFmt = $self->getProperty('InputAnnotFormat');
    my $inputContig = $self->getProperty('InputContigFile');
    my $inputContigFmt = $self->getProperty('InputContigFormat');
    my $transcriptTag = $self->getProperty('transcriptTag');
    my $cdsTag = $self->getProperty('CDSTag');
#    my $coverageDensity = $self->getProperty('coverageDensity');
    my $gcClasses = $self->getProperty('gcClasses');
    my $geneModelType = $self->getProperty('geneModelType');
    my $numDJobNodes = $self->getProperty('numDJobNodes');
    my $djobNodeClass = $self->getProperty('DJobNodeClass'); 
    my $prefixOutputDir = $self->getProperty('prefixOutputDir');
    my $DJobInputBaseDir = $prefixOutputDir.".test";
    my $species = $self->getProperty('species');
    my $inputAnnotFile = $self->getProperty('InputAnnotFile');
    my $inputContigFile = $self->getProperty('InputContigFile');
    
    `mkdir $DJobInputBaseDir` unless -d $DJobInputBaseDir;
    `mkdir $DJobInputBaseDir/preprocess` unless -d "$DJobInputBaseDir/preprocess";
    `rm -rf $DJobInputBaseDir/preprocess/*`; 
    `mkdir $DJobInputBaseDir/model` unless -d "$DJobInputBaseDir/model";
    `rm -rf $DJobInputBaseDir/model/*`; 

    my $modelOutputDir = $prefixOutputDir.".model";
    my $minPercJunctionAlignments = $self->getProperty('minPercJunctionAlignments');
    
    my $genUTROnlyModel = $self->getProperty('genUTROnlyModel');
    my $modelUTRs = $self->getProperty('modelUTRs');
    my $trainingIterations = $self->getProperty('trainingIterations');
    my $closeSpecies = $self->getProperty("closestSpecies");
    my $cmd = "$ENV{CRAIG_HOME}/python/bin/craigReannotate.py --output-model-scores --pre-config $inputDir/preconfig --prefix-output-dir $prefixOutputDir --annot-fmt $inputAnnotFmt --contig-fmt $inputContigFmt --transcript-tag $transcriptTag --cds-tag $cdsTag --gc-classes $gcClasses --training-iterations $trainingIterations --force-train --model $geneModelType --min-perc-junc-aligns $minPercJunctionAlignments --djob-num-nodes  $numDJobNodes --djob-node-class $djobNodeClass --djob-input $DJobInputBaseDir";
    
    $cmd = $cmd." --closest-species $closeSpecies" if $closeSpecies ne "none";
    $cmd = $cmd." --utr-only-model" if $genUTROnlyModel == "yes";
    $cmd = $cmd." --model-utrs" if$modelUTRs == "yes";
    $cmd = $cmd." $species $inputAnnotFile $inputContigFile";
#    exit(0);  ##for testin
    print "CRAIG Run Command: $cmd\n";
    $self->{nodeForInit}->runCmd($cmd);

}

sub initNode {
    my ($self, $node, $inputDir) = @_;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;
    return 1;
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;
    my $cmd = "echo nop";
    return $cmd;
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
}

sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;
  return 1;
}

1;
