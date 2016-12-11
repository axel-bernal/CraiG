package DJob::DistribJobTasks::XSignalsTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["prefixFilesPath", "",  "prefix to RNA-seq junction and coverage files"],
 ["inputFilePath", "",  "full path to input FASTA file"],
 ["blockLen", "20",  "block length used in discarding permutations [20]"],
 ["chunkSize", "30000",  "fragment size used to partition contigs larger than chunkSize [30000]"],
 ["minCoverage", "0", "Minimum coverage required. If given, regions below that threshold will be ignored. This option speeds up the program dramatically and it is highly recommended to be set (typically to a value of 2) when working on entire chromosomes/long contigs [0]"],
 ["stranded", "no", "whether input RNA-Seq data is stranded or not yes [no]\n"],
 ["numPermutations", "1000", "number of permutations for computing change points [1000]"]
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

sub initServer {
    my ($self, $inputDir) = @_;
    my $prefixOutputFiles = $self->getProperty("prefixFilesPath");
    my $fastaFileName = $self->getProperty("inputFilePath");
    my $chunkSize = $self->getProperty("chunkSize");
    
    my $date = `date`;
    open(LOGFILE, ">$inputDir/xsignals.log_master") || die "Can't open file $inputDir/xsignals.log_master for writing\n";
    print LOGFILE "\nstart: $date\n";
    print LOGFILE "Splitting contigs\n";
    my $cmd = "$ENV{CRAIG_HOME}/perl/bin/splitContigs.pl $chunkSize < $fastaFileName > $inputDir/split_contigs.fa";
    print LOGFILE $cmd, "\n";

#    $self->{nodeForInit}->runCmd("rm -f $prefixOutputFiles.xsignal");
#    $self->{nodeForInit}->runCmd("rm -f $prefixOutputFiles.xsigscores");
    $self->{nodeForInit}->runCmd("rm -f $prefixOutputFiles.stderr");
    close LOGFILE;
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

}

sub getInputSetSize {
    my ($self, $inputDir) = @_;
    my $fastaFileName = $self->getProperty("inputFilePath");
    my $chunkSize = $self->getProperty("chunkSize");

    $self->{nodeForInit}->runCmd("gunzip $fastaFileName.gz") unless !(-e "$fastaFileName.gz");

    my $cmd = "$ENV{CRAIG_HOME}/perl/bin/splitContigs.pl $chunkSize < $fastaFileName > $inputDir/split_contigs.fa";
    $self->{nodeForInit}->runCmd($cmd);
    print "Counting sequences in $inputDir/split_contigs.fa\n";
    $self->{fastaFile} = CBIL::Bio::FastaFileSequential->new("$inputDir/split_contigs.fa");
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir,$subTask) = @_;

    if(!$subTask->getRedoSubtask()){
	$self->{fastaFile}->writeSeqsToFile($start, $end, 
					    "$subTaskDir/seqsubset.fa");
    }

    $node->runCmd("touch $subTaskDir/seqsubset.fa.touch",1);
    $node->runCmd("/bin/rm $subTaskDir/seqsubset.fa.touch",1);
    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $blockLen = $self->getProperty("blockLen");
    my $prefixFilesPath = $self->getProperty("prefixFilesPath");
    my $stranded = ($self->getProperty("stranded") eq "no"?"":"--stranded");
    my $chunkSize = $self->getProperty("chunkSize");
    my $num_perms = $self->getProperty("numPermutations");
    my $minCoverage = $self->getProperty("minCoverage");
    my $cmd = "$ENV{CRAIG_HOME}/bin/get_covxsignals $stranded --min-cov=$minCoverage --block-len=$blockLen --chunk-size=$chunkSize --prefix-evidence=$prefixFilesPath --prefix-output=seqsubset --num-permutations=$num_perms seqsubset.fa";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    my $prefixOutputFiles = $self->getProperty("prefixFilesPath");
    $node->runCmd("cat $nodeExecDir/seqsubset.xsignal | $ENV{CRAIG_HOME}/perl/bin/fastalike2files.pl --append -prefix $mainResultDir/split.xsignal");
    $node->runCmd("cat $nodeExecDir/seqsubset.xsigscores | $ENV{CRAIG_HOME}/perl/bin/fastalike2files.pl --append -prefix $mainResultDir/split.xsigscores");
    return 1 if $node->getErr();
    $node->runCmd("cat $nodeExecDir/subtask.stderr >> $prefixOutputFiles.stderr");
}

sub cleanUpServer {
    my($self, $inputDir, $mainResultDir, $node) = @_;
    my $prefixOutputFiles = $self->getProperty("prefixFilesPath");
    print "Postprocessing output\n";
    my $cmd = "$ENV{CRAIG_HOME}/perl/bin/postProcessXSignalTask.pl --mainResultDir $mainResultDir --output output --numContigs ".$self->{fastaFile}->getCount();
    $node->runCmd($cmd);
    print "$cmd is done!\n";
    my $buff_size = $self->{fastaFile}->getCount() + 1000;
    print "Removing temporary files .. ";
    $cmd = "find $mainResultDir -iname \"split.*\" -print0 | xargs -0 -s $buff_size rm -f";
    $node->runCmd($cmd);
    print "done!\n";
    $cmd  = "cp $mainResultDir/output.xsignal $prefixOutputFiles.xsignal";
    $node->runCmd($cmd);
    $cmd  = "cp $mainResultDir/output.xsigscores $prefixOutputFiles.xsigscores";
    $node->runCmd($cmd);

}

1;
