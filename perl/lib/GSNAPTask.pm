package DJob::DistribJobTasks::GSNAPTask;
use DJob::DistribJob::Task;
use File::Basename;
use File::Path;
use Cwd;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["readFilePath",   "",     "full path to read file"],
 ["pairedReadFilePath",   "none",     "full path to paired read file (optional)"],
 ["genomeDatabase",   "",     "genome database id"],
 ["genomeFastaFile",   "",     "genome file in fasta format"],
 ["geneAnnotationFile",   "none",     "geneAnnotationFile in gff3 format"],
 ["CDSTag", "CDS", "Tag for recognizing coding regions in InputAnnotFile (exon | [CDS])"],
 ["createSAMFile",   "false",     "create SAM file if true ([false] | true)"],
 ["orientation",     "FR",    "orientation of paired or single reads([FR]| RF | FF | NN | R | F | N)"],
 ["createJunctionsFile", "false", "create Junctions file if true ([false] | true)"],
 ["createBigWigFiles", "false", "create Coverage BigWig files if true ([false] | true)"],
 ["SNPs",   "false",     "run snp finder -- not implemented yet -- ([false] | true)"],
 ["postProcess",   "false",     "run Post Processing steps if true ([false] | true)"],
 );

## note passing in additional param that skips computing the size in constructor
sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties, 1); 
    return $self;
}

# called once .. do everything here upstream of sending chunks to nodes
# if have pairedReadFile check to see that size of chunk is even number
# combine files into single fasta file
# get qual values?
# either create fasta index here or in getCount method
sub initServer {
    my ($self, $inputDir) = @_;
    
    ## NOTE: need to make certain that subtask size is an even number.
    die "subTaskSize must be an even number!\n" if ($self->{subTaskSize} % 2);
    
    my $readFilePath = $self->getProperty("readFilePath");
    my $genomeFastaFile = $self->getProperty("genomeFastaFile");
    my $geneAnnotationFile = $self->getProperty("geneAnnotationFile");
    my $pairedReadFilePath = $self->getProperty("pairedReadFilePath");
    my $genomeDB = $self->getProperty("genomeDatabase");
    my $createSAMFile = $self->getProperty("createSAMFile");
    my $cdsTag = $self->getProperty("CDSTag");

    die "readFilePath $readFilePath does not exist" unless -e "$readFilePath";
    die "pairedReadFilePath $pairedReadFilePath does not exist" if $pairedReadFilePath ne 'none' && !(-e "$pairedReadFilePath");
    die "readFilePath equals pairedReadFilePath" if $pairedReadFilePath != 'none' && $pairedReadFilePath eq $readFilePath;
    die "genomeFastaFile $genomeFastaFile does not exist" unless -e "$genomeFastaFile";

    ##if aren't creating SAM file then don't need to create qual file
    if($createSAMFile =~ /true/i){
	$self->{quals} = 1;
    }else{
	$self->{quals} = 0;
    }

    # creating length file
    open(IC, "<$genomeFastaFile") || die "Cannot read from $genomeFastaFile\n";
    open(OC, ">$inputDir/contigs_tmp.fa") || die "Cannot write in $inputDir/contigs_tmp.fa\n";
    my $line;
    my $numContigs = 0;
    while(defined($line = <IC>)) {
	chomp($line);
	if($line =~ /^>(\S+)\s*/) {  print OC "\n>$1\n"; $numContigs++; } 
	else { print OC $line; }
    }
    close(OC);
    $self->{numContigs} = $numContigs;
    $self->{nodeForInit}->runCmd("cat $inputDir/contigs_tmp.fa | sed 1,1d | sed -e '\$a\\' > $inputDir/contigs.fa && rm -f $inputDir/contigs_tmp.fa");
    $self->{nodeForInit}->runCmd("$ENV{CRAIG_HOME}/perl/bin/contigLength.pl -ctg $inputDir/contigs.fa > $inputDir/contigs.len");
    ## initializing database or genome directory if fasta and annotation files are provided
    my $dirDatabases = $self->{nodeForInit}->runCmd("gmap -d \'?\' | grep Available | cut -f6 -d \" \" | cut -f1 -d \":\"");
    chomp($dirDatabases);

    if(!(-e "$dirDatabases/$genomeDB/$genomeDB.maps/$genomeDB.splicesites.iit")) {
	print "Generating GSNAP indexes for genome and gene annotation...\n";
	die "geneAnnotationFile $geneAnnotationFile does not exist" unless -e "$geneAnnotationFile";
	rmtree("$inputDir/contigs") unless !(-d "$inputDir/contigs");
	rmtree("$dirDatabases/$genomeDB") unless !(-d "$dirDatabases/$genomeDB");
	$self->{nodeForInit}->runCmd("$ENV{CRAIG_HOME}/perl/bin/formatFastaFile.pl < $inputDir/contigs.fa | $ENV{CRAIG_HOME}/perl/bin/fastalike2files.pl -dir $inputDir/contigs");
	$self->{nodeForInit}->runCmd("ls $inputDir/contigs/* > $inputDir/contigs.list");
	$self->{nodeForInit}->runCmd("gmap_build -d $genomeDB -k 15 $inputDir/contigs.list");
	$self->{nodeForInit}->runCmd("cat $geneAnnotationFile | grep -P \"\t$cdsTag\t\" | $ENV{CRAIG_HOME}/perl/bin/gff32GTF.pl -all | sed 's/\tCDS\t/\texon\t/g' | gtf_splicesites > $inputDir/$genomeDB.splicesites");
	$self->{nodeForInit}->runCmd("cat $inputDir/$genomeDB.splicesites | iit_store -o $inputDir/$genomeDB.splicesites.iit");
	$self->{nodeForInit}->runCmd("cp $inputDir/$genomeDB.splicesites.iit $dirDatabases/$genomeDB/$genomeDB.maps");
	print " Done!\n";
    }

    die "genome database $dirDatabases/$genomeDB does not exist" unless -e "$dirDatabases/$genomeDB";

    ##make fasta file and qual files .. note that will chunk them up here as well.
    my $date;
    if(!(-d "$inputDir/subtasks")){
	mkdir("$inputDir/subtasks");
    }
    my $madeReadFile = 0;
  REDOSUBTASKS:
    if(!(-e "$inputDir/subtasks/sequenceCount1") || 
       ((-e $pairedReadFilePath) && !(-e "$inputDir/subtasks/sequenceCount2"))) {  ##file will exist if have already done this
      $date = `date`; chomp $date;
      print "[$date] parsing reads file(s) to fasta and qual files\n";
#      $self->{nodeForInit}->runCmd("cat $inputDir/quals | makeSubtasksFromFAFileForRUM.pl --stdin --outStem quals --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize}") if $self->{quals};
      $self->{nodeForInit}->runCmd("$ENV{CRAIG_HOME}/perl/bin/makeSubtasksFromFastQ.pl --fileName $readFilePath --outStem reads1 --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize} > $inputDir/subtasks/sequenceCount1");
      if(-e $pairedReadFilePath) {
	  $self->{nodeForInit}->runCmd("$ENV{CRAIG_HOME}/perl/bin/makeSubtasksFromFastQ.pl --fileName $pairedReadFilePath --outStem reads2 --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize} > $inputDir/subtasks/sequenceCount2");
      }
      $madeReadFile = 1;
    }

    my $seqCount1 = `cat $inputDir/subtasks/sequenceCount1`;
    if(-e $pairedReadFilePath) { 
	my $seqCount2 = `cat $inputDir/subtasks/sequenceCount2`;
	if($seqCount1 ne $seqCount2) {
	    print "Error: Paired Input set sizes are different: $seqCount1 $seqCount2\n";
	    exit 1;
	}
    }
    ##set the count .. inputSetSize
    chomp $seqCount1;
    $self->{inputSetSize} = $seqCount1;

    ##get and report size ....
    $self->{size} = $self->getInputSetSize($inputDir);
    if ($self->{size} == 0) {
	print "Error: Input set size is 0\n";
	exit 1;
    }
    my $c = int($self->{size} / $self->{subTaskSize});
    $c += 1 if $self->{size} % $self->{subTaskSize};
    $self->{numSubtasks} = $c;

    ##check to see if the the number of files created is consistent with number of subtasks
    ## if not, make the files above because subtask size has changed;
    ## not going to test for change in inputfile as if user is silly enough to change insitu they'll have to live with it
    ## intended to catch case where subtask size has been changed.
    my @sFiles = glob("$inputDir/subtasks/reads1.*");
    if(scalar(@sFiles) != $self->{numSubtasks}){
      die "Unable to create subtask reads files properly\n" if $madeReadFile;
      system("/bin/rm $inputDir/subtasks/*");
      goto REDOSUBTASKS;
    }
    print "Input set size is $self->{size}\n";
    print "Subtask set size is $self->{subTaskSize} ($self->{numSubtasks} subtasks)\n";
    
    $self->{pairedEnd} = (-e "$pairedReadFilePath") ? "paired" : "single";
               
    $date = `date`;
    open(LOGFILE, ">gsnap.log_master");
    print LOGFILE "\nstart: $date\n";
    print LOGFILE "paired- or single-end: $self->{pairedEnd}\n";
    print LOGFILE "num subtasks: $self->{numSubtasks}\n";
    close(LOGFILE);
    
}

sub initNode {
  my ($self, $node, $inputDir) = @_;
  return 1;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    return $self->{inputSetSize};
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;

    my $subtaskNum = int($start / $self->{subTaskSize}) + 1;
    $node->runCmd("touch $inputDir/subtasks/reads1.$subtaskNum.touch",1);
    $node->runCmd("/bin/rm $inputDir/subtasks/reads1.$subtaskNum.touch",1);
    $node->runCmd("cp $inputDir/subtasks/reads1.$subtaskNum $nodeExecDir/seqSubset1.fa");
    $node->runCmd("cp $inputDir/contigs.len $nodeExecDir");

    if($self->{pairedEnd} eq "paired") {
	$node->runCmd("cp $inputDir/subtasks/reads2.$subtaskNum $nodeExecDir/seqSubset2.fa");
    }	
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir,$subtaskNumber,$mainResultDir) = @_;
    
    my $genomeDB = $self->getProperty("genomeDatabase");
    my $createSAMFile = $self->getProperty("createSAMFile");
    my $orientation = $self->getProperty("orientation");
    my $createJunctionsFile = $self->getProperty("createJunctionsFile");
    my $createBigWigFiles = $self->getProperty("createBigWigFiles");

    my $cmd = "$ENV{CRAIG_HOME}/perl/bin/runGSNAPOnNode.pl --genomeDB $genomeDB --readseq seqSubset1.fa".($self->{pairedEnd} eq "paired" ? " --paired-readseq seqSubset2.fa " : "").($createJunctionsFile eq "true" ? " --createJunctionsFile " : "").($createBigWigFiles eq "true" ? " --createBigWigFiles " : "").($createSAMFile eq "true" ? " --createSAMFile " : "")." --contig-len contigs.len --orientation $orientation --subTaskNumber $subtaskNumber --output $mainResultDir/GSNAP";
# -o $orientation .. not using this for the moment
    
    $self->{subtaskCmd} = $cmd;
    return $self->{subtaskCmd};
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    return 1 if $node->getErr();

}


# do all post-processing here:
# a node is now passed in that can be used to run commands on a node using $node->runCmd("cmd")
# NOTE that in order to use this must add keepNodeForPostProcessing=yes to controller.prop file

sub cleanUpServer {
    my($self, $inputDir, $mainResultDir, $node) = @_;

    my $genomeDB = $self->getProperty("genomeDatabase");
    my $createSAMFile = $self->getProperty("createSAMFile");
    my $orientation = $self->getProperty("orientation");
    my $createBigWigFiles = $self->getProperty("createBigWigFiles");
    my $createJunctionsFile = $self->getProperty("createJunctionsFile");
    my $stNum;

    unlink("$mainResultDir/GSNAP.aln");
    for($stNum = 1; $stNum <= $self->{numSubtasks}; $stNum++) {
	my $outfile = "$mainResultDir/GSNAP.$stNum.aln";
	die "file $outfile cannot be found" unless -e "$outfile";
	$node->runCmd("cat $outfile >> $mainResultDir/GSNAP.aln");
    }

    if($createSAMFile eq "true") {
	unlink("$mainResultDir/GSNAP.sam");
	for($stNum = 1; $stNum <= $self->{numSubtasks}; $stNum++) {
	    my $outfile = "$mainResultDir/GSNAP.$stNum.sam";
	    die "file $outfile cannot be found" unless -e "$outfile";
	    $node->runCmd("cat $outfile >> $mainResultDir/GSNAP.sam");
	}
    }
	
    
    return 1 if $self->getProperty('postProcess') ne "yes";

    print "Postprocessing output for $genomeDB\n";
    my $cmd = "$ENV{CRAIG_HOME}/perl/bin/postProcessGSNAPTask.pl".($createJunctionsFile eq "true" ? " --createJunctionsFile " : "").($createBigWigFiles eq "true" ? " --createBigWigFiles " : "")." --mainResultDir $mainResultDir --orientation $orientation --output GSNAP --numContigs $self->{numContigs} --numSubTasks $self->{numSubtasks}";
    $node->runCmd($cmd);
    print "$cmd is done!\n";
    my $buff_size = $self->{numSubtasks} + 1000;
    print "Removing temporary files .. ";
    $cmd = "find $mainResultDir -iname \"GSNAP.*.al*\" -print0 | xargs -0 -s $buff_size rm -f";
    $node->runCmd($cmd);
    $cmd = "find $mainResultDir -iname \"GSNAP.*.Unique.cov\" -print0 | xargs -0 -s $buff_size rm -f";
    $node->runCmd($cmd);
    if($createSAMFile eq "true") {
	$cmd = "find $mainResultDir -iname \"GSNAP.*.sam\" -print0 | xargs -0 -s $buff_size rm -f";
	$node->runCmd($cmd);
    }

    $buff_size = $self->{numContigs} + 1000;
    $cmd = "find $mainResultDir -iname \"GSNAP.Unique.cov.*\" -print0 | xargs -0 -s $buff_size rm -f";
    $node->runCmd($cmd);
    print "Postprocessing Done!\n";
    return 1;
}

1;
