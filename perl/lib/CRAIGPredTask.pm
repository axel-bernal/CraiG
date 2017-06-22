package DJob::DistribJobTasks::CRAIGPredTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["paramsFile", "", "Name of the file containing the gene model parameters"],
 ["prefixFilesPath", "",  "prefix to resource files, such as RNA-seq junction, coverage files and/or evidence signal and state files for ensemble models"],
 ["inputFilePath", "",  "full path to input FASTA file"],
 ["predictUTRs", "yes",  "Predictions will contain UTR regions. If not specified and paramters model UTRs, expect worse performance, so choose carefully"],
 ["format", "locs", "output format"],
 ["START-resource", " ", "resource name for Start signals"],
 ["STOP-resource", " ", "resource name for Stop signals"],
 ["DONOR-resource", " ", "resource name for Donor signals"],
 ["ACCEPTOR-resource", " ", "resource name for Acceptor signals"]
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

sub initServer {
    my ($self, $inputDir) = @_;
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->getProperty("inputFilePath");

    if (-e "$fastaFileName.gz") {
	&runCmd("gunzip $fastaFileName.gz");
    }

    print "Counting sequences in $fastaFileName\n";
    $self->{fastaFile} = CBIL::Bio::FastaFileSequential->new($fastaFileName);
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir,$subTask) = @_;

    if(!$subTask->getRedoSubtask()) {
      $self->{fastaFile}->writeSeqsToFile($start, $end, 
					  "$subTaskDir/seqsubset.fa");
    }

    $node->runCmd("touch $subTaskDir/seqsubset.fa.touch",1);
    $node->runCmd("/bin/rm $subTaskDir/seqsubset.fa.touch",1);
    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $format = $self->getProperty("format");
    my $paramsFile = $self->getProperty("paramsFile");
    my $predictUTRs = ($self->getProperty("predictUTRs") eq "yes" ? "--predict-utrs" : "");
    my $prefixFilesPath = $self->getProperty("prefixFilesPath");
    my $startResource = ($self->getProperty("START-resource") eq " " ? "" : "--START-resource=".$self->getProperty("START-resource"));
    my $stopResource = ($self->getProperty("STOP-resource") eq " " ? "" : "--STOP-resource=".$self->getProperty("STOP-resource"));
    my $donorResource = ($self->getProperty("DONOR-resource") eq " " ? "" : "--DONOR-resource=".$self->getProperty("DONOR-resource"));
    my $acceptorResource = ($self->getProperty("ACCEPTOR-resource") eq " " ? "" : "--ACCEPTOR-resource=".$self->getProperty("ACCEPTOR-resource"));
    my $cmd = "$ENV{CRAIG_HOME}/bin/craigPredict $predictUTRs  $startResource $stopResource $donorResource $acceptorResource --output-file=seqsubset.locs --format=$format --best=1 --strand=both --prefix-evidence=$prefixFilesPath $paramsFile seqsubset.fa";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    $node->runCmd("cat $nodeExecDir/seqsubset.locs > $mainResultDir/seqsubset$subTaskNum.locs");
    $node->runCmd("cat $nodeExecDir/subtask.stderr >> $mainResultDir/stderr0");
    return 1 if $node->getErr();

}


sub cleanUpServer {
    my($self, $inputDir, $mainResultDir, $node) = @_;
    my $gen_count = 0;
    my $tr_count = 0;
    my $numSubTasks = 1;
    open(ALLPRED, ">>$mainResultDir/output0.locs");
    while(1) {
	last unless open(NODEPRED, "<$mainResultDir/seqsubset$numSubTasks.locs");
	my $line = <NODEPRED>;
	my $currgen_count = $gen_count;
	my $currtr_count = $tr_count;
	while(defined($line) && $line =~ /^>gen_(\d+)\|tr_(\d+)\s+(\S+)$/) {
	    $gen_count = $currgen_count + $1;
	    $tr_count = $currtr_count + $2;
	    print ALLPRED ">gen_".$gen_count."|tr_".$tr_count."\t$3\n";
	    $line = <NODEPRED>;
	}

#	unlink("$mainResultDir/seqsubset$numSubTasks.locs");
	$numSubTasks++;
	close(NODEPRED);
    }
    close(ALLPRED);
}

1;
