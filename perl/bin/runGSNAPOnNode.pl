#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use File::Temp;
use strict;
use Getopt::Long;

my $genomeDB = "";
my $readseq = "";
my $paired_readseq = "";
my $output = "";
my $createSAMFile = 0;
my $createJunctionsFile = 0;
my $createBigWigFiles = 0;
my $subTaskNumber = 0;
my $orientation = "";
my $contiglen_file = "";

&GetOptions("-genomeDB=s" => \$genomeDB,
	    "-readseq=s" => \$readseq,
	    "-paired-readseq=s" => \$paired_readseq,
	    "-output=s" => \$output,
	    "-createSAMFile" => \$createSAMFile,
	    "-createJunctionsFile" => \$createJunctionsFile,
	    "-createBigWigFiles" => \$createBigWigFiles,
	    "-subTaskNumber=i" => \$subTaskNumber,
	    "-orientation=s" => \$orientation,
	    "-contig-len-file=s" => \$contiglen_file,
    );

die "Must provide required parameters" unless $genomeDB ne "" && $output ne "" && $readseq ne "" && $subTaskNumber > 0 && $orientation ne "" && $contiglen_file ne "";

my (undef, $filename1) = File::Temp::tempfile();
print "command to run: gsnap -d $genomeDB -k 15 -s $genomeDB.splicesites.iit -N 1 -t 2 $readseq $paired_readseq > $filename1\n";
`gsnap -d $genomeDB -k 15 -s $genomeDB.splicesites.iit -N 1 -t 2 $readseq $paired_readseq > $filename1`;
`cp $filename1 $output.$subTaskNumber.aln`;

if($createJunctionsFile) {
    my (undef, $filename2) = File::Temp::tempfile();
    `cat \"$output.$subTaskNumber.aln\" | junctions2weightedLocs.py --dataset-id $subTaskNumber --orientation $orientation --format gsnap > $filename2`;
    `cp $filename2 $output.$subTaskNumber.alljunction.locs`;
}

if($createBigWigFiles) {
    my (undef, $filename3) = File::Temp::tempfile();
    `cat \"$output.$subTaskNumber.aln\" | gsreads2covfilter.py --orientation $orientation --len-file $contiglen_file > $filename3`;
    `cp $filename3 $output.$subTaskNumber.Unique.cov`;
}

if($createSAMFile) {
    my (undef, $filename4) = File::Temp::tempfile();
    print "command to run: gsnap -d $genomeDB -k 15 -s $genomeDB.splicesites.iit -N 1 -t 2 $readseq $paired_readseq > $filename4\n";
    `gsnap --no-sam-headers -A sam -d $genomeDB -k 15 -s $genomeDB.splicesites.iit -N 1 -t 2 $readseq $paired_readseq > $filename4`;
    `cp $filename4 $output.$subTaskNumber.sam`;
}
