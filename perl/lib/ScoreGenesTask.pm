package DJob::DistribJobTasks::ScoreGenesTask;

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
 ["prefixEvidencePath", "",  "prefix path to RNA-seq junction and coverage files"],
 ["fastaFilePath", "",  "full path to input FASTA file"],
 ["annotFilePath", "", "Full path to genome isoform annotation"],
 ["haveUTRs", "", "Set to 'yes' if genes have UTRs otherwise 'no' [yes]"]
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

sub initServer {
    my ($self, $inputDir) = @_;
    my $date = `date`;
    open(LOGFILE, ">$inputDir/scoregenes.log_master") || die "Can't open file $inputDir/scoregenes.log_master for writing\n";
    print LOGFILE "\nstart: $date\n";
    print LOGFILE "Creating subcontigs\n";
    close LOGFILE;
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

}

sub getInputSetSize {
    my ($self, $inputDir) = @_;
    my $fastaFilePath = $self->getProperty("fastaFilePath");
    my $annotFilePath = $self->getProperty("annotFilePath");

    $self->{nodeForInit}->runCmd("gunzip $fastaFilePath.gz") unless !(-e "$fastaFilePath.gz");

    my $cmd = "cat $annotFilePath | $ENV{CRAIG_HOME}/perl/bin/extractIGenics.pl -ctg $fastaFilePath > $inputDir/input_subcontigs.igenics";    
    $self->{nodeForInit}->runCmd($cmd);
    $cmd = "$ENV{CRAIG_HOME}/perl/bin/getSubContigsFromIgenics.pl --use-subcontiglocs-as-ids < $inputDir/input_subcontigs.igenics | $ENV{CRAIG_HOME}/perl/bin/filterRepeatedGenes.pl > $inputDir/input_subcontigs.subcontigs";
    $self->{nodeForInit}->runCmd($cmd);
    $cmd = "$ENV{CRAIG_HOME}/perl/bin/moveLocsToSubContigs.pl --preserve-ids -l -ctg $inputDir/input_subcontigs.subcontigs < $annotFilePath > $inputDir/input_subcontigs.locs";
    $self->{nodeForInit}->runCmd($cmd);
    $cmd = "$ENV{CRAIG_HOME}/perl/bin/getSeqsFromGeneLocs.pl -ctg $fastaFilePath < $inputDir/input_subcontigs.subcontigs | cut -f1 > $inputDir/input_subcontigs.fa";
    $self->{nodeForInit}->runCmd($cmd);
    print "Counting sequences in $inputDir/input_subcontigs.fa\n";
    $self->{fastaFile} = CBIL::Bio::FastaFileSequential->new("$inputDir/input_subcontigs.fa");
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir, $subTask) = @_;
    
    if(!$subTask->getRedoSubtask()) {
	$self->{fastaFile}->writeSeqsToFile($start, $end, 
					    "$subTaskDir/seqsubset.fa");
    }
    
    &runCmd("touch $subTaskDir/seqsubset.fa.touch",1);
    &runCmd("/bin/rm $subTaskDir/seqsubset.fa.touch",1);
    &runCmd("$ENV{CRAIG_HOME}/python/bin/biointers.py -by ci $inputDir/input_subcontigs.locs $subTaskDir/seqsubset.fa > $subTaskDir/seqsubset.locs");
    $node->runCmd("touch $subTaskDir/seqsubset.locs.touch",1);
    $node->runCmd("/bin/rm $subTaskDir/seqsubset.locs.touch",1);    
    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;
    my $paramsFile = $self->getProperty("paramsFile");
    my $haveUTRs = $self->getProperty("haveUTRs");
    my $prefixEvidencePath = $self->getProperty("prefixEvidencePath");
    $node->runCmd("sync");
    my $cmd = "$ENV{CRAIG_HOME}/bin/score_genes --output=seqsubset.gms.locs ".($haveUTRs eq "yes" ? "--have-utrs ": "")."--prefix-evidence=$prefixEvidencePath $paramsFile seqsubset.fa seqsubset.locs";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    $node->runCmd("cat $nodeExecDir/seqsubset.gms.locs >> $mainResultDir/input_subcontigs.gms.locs");
    return 1 if $node->getErr();
    $node->runCmd("cat $nodeExecDir/subtask.stderr >> $mainResultDir/input_subcontigs.gms.stderr");
}

sub cleanUpServer {
    my($self, $inputDir, $mainResultDir, $node) = @_;
    print "Postprocessing output\n";
    my $cmd = "$ENV{CRAIG_HOME}/perl/bin/moveLocsToOrigContigs.pl $inputDir/input_subcontigs.subcontigs < $mainResultDir/input_subcontigs.gms.locs > $mainResultDir/output.gms.locs";
    $node->runCmd($cmd);
    print "$cmd is done!\n";
}

1;
