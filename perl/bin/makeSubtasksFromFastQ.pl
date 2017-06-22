#!/usr/bin/perl

use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use Getopt::Long;

my $debug = 0;

$| = 1;
my $firstNArow = 0;
my $secondNArow = 0;
my ($fileName,$directory,$subtaskSize);
my $outStem = "subTask.fa";

&GetOptions("fileName=s" => \$fileName, ##just one fasta file as input
            "outStem=s" => \$outStem,
            "subtaskSize=i" => \$subtaskSize,        ##take fasta on stdin
            "directory=s" => \$directory  ##directory for output
            );

die "invalid directory '$directory'" unless -d "$directory";
die "you must specify --fileName <inputfile>" unless (-e $fileName);
die "you must specify --subtaskSize <size>" unless $subtaskSize;

my $standard = "true";
open(INFILE, $fileName);
for(my $i=0; $i<1000; $i++) {
    my $line = <INFILE>;
    chomp($line);
    if($line eq '') {
        if($i < 4) {
            $standard = "false";
        }
	$i = 1000;
    } else  {
	if(!($line =~ /^@/)) {
	    $standard = "false";
	}
	$line = <INFILE>;
	chomp($line);
	if(!($line =~ /^[ACGTN.]+$/)) {
	    $standard = "false";
	}
	$line = <INFILE>;
	chomp($line);
	if(!($line =~ /^\+/)) {
	    $standard = "false";
	}
	$line = <INFILE>;
	chomp($line);
	if(!($line =~ /^[!-~]+$/)) {
	    $standard = "false";
	}
    }
}
close(INFILE);

if($standard eq 'true') {
    $firstNArow = 0;
    $secondNArow = 4;
} else {
    print "\nSorry, can't figure these files out, doesn't look like fastq...\n";
    exit(0);
}

open(F,"$fileName") || die "$fileName cannot be opened\n";
writeStems(\*F, $directory, $outStem, $firstNArow, $secondNArow, $subtaskSize);
close(F);

sub writeStems {
    my($in_fh, $directory, $outStem, $firstNArow, $secondNArow, $subtaskSize) = @_;
    my $ct = 0;
    my $cnt = 0;
    my $linecnt = 0;
    my $stnum = 1;
    my $block = $secondNArow - $firstNArow;
    open(OUT,">$directory/$outStem.$stnum");
    my $line;
    while(defined($line = <$in_fh>)) {
	$linecnt++;
	if((($cnt - $firstNArow) % $block) == 0) {
	    if($line !~ /^@/) {
		print STDERR "\nWARNING: in script makeSubtasksFronFastQ.pl: There seems to be something wrong with line $linecnt in input file\nIt should be an id line  but instead it is:\n$line\n\n";
		print "\nWARNING: in script makeSubtasksFronFastQ.pl: There seems to be something wrong with line $linecnt in input file\nIt should be an id line  but instead it is:\n$line\n\n";
		print "\nSorry, can't figure these files out, maybe they're not fastq, they're corrupt, or you might have to write a custom parser.\n";
		exit();
	    }
	    $ct++;
	    if($ct != 1 && $ct % $subtaskSize == 1){
		close(OUT);
		$stnum++;
		open(OUT,">$directory/$outStem.$stnum");
	    }
	}
	$cnt++;
	print OUT $line;
    }
    
    if($linecnt % $block != 0) {
	print STDERR "\nWarning: in script makeSubtasksFromFastQ.pl: the last block of lines in input fastq file is not the right size.\n\n";
    }    
##print total number of sequences
    print "$ct\n";
}


