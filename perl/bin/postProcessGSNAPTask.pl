#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use File::Temp;
use strict;
use Getopt::Long;

my $output = "";
my $createSAMFile = 0;
my $createJunctionsFile = 0;
my $createBigWigFiles = 0;
my $numContigs = 0;
my $numSubTasks = 0;
my $mainResultDir = "";
my $orientation = "N";

&GetOptions("-output=s" => \$output,
	    "-mainResultDir=s" => \$mainResultDir,
	    "-createSAMFile" => \$createSAMFile,
	    "-createJunctionsFile" => \$createJunctionsFile,
	    "-createBigWigFiles" => \$createBigWigFiles,
	    "-numContigs=i" => \$numContigs,
	    "-numSubTasks=i" => \$numSubTasks,
	    "-orientation=s" => \$orientation,
    );

die "Must provide required parameters" unless  $output ne "" &&  $numContigs > 0 && $mainResultDir ne "" && $numSubTasks > 0;

my (undef, $filename1) = File::Temp::tempfile();
my (undef, $filename2) = File::Temp::tempfile();
my $stNum;

if($createJunctionsFile) {
    print "Creating junctions $mainResultDir/$output.junction.locs .. "; 
    if(!(-e "$mainResultDir/$output.junction.locs")) {
	for($stNum = 1; $stNum <= $numSubTasks; $stNum++) {
	    my $outfile = "$mainResultDir/$output.$stNum.alljunction.locs";
	    die "file $outfile cannot be found" unless -e "$outfile";
	    `cat \"$outfile\" >> $filename2`;
	}
	`cp $filename2 $mainResultDir/$output.alljunction.locs`;
	`filterRepeatedGenes.pl < $mainResultDir/$output.alljunction.locs > $filename1`;
	`cp $filename1 $mainResultDir/$output.junction.locs`;
    }
    print " done!\n";
}
    
my (undef, $filename3) = File::Temp::tempfile();
my $buff_size = $numContigs + 1000;

if($createBigWigFiles) {
    print "Creating coverage $mainResultDir/$output.Unique.cov .. ";
    if(!(-e "$mainResultDir/$output.Unique.cov")) {
	`find $mainResultDir -iname \"$output.Unique.cov.*\" -print0 | xargs -0 -s $buff_size rm -f`;
	for($stNum = 1; $stNum <= $numSubTasks; $stNum++) {
	    `cat \"$mainResultDir/$output.$stNum.Unique.cov\" | fastalike2files.pl --append -prefix $mainResultDir/$output.Unique.cov`;
	}

	my @outfiles = `find $mainResultDir -iname \"$output.Unique.cov.*\" -print0 | xargs -0 -s $buff_size ls`;
	foreach my $outfile (@outfiles) {
	    chomp $outfile;
	    die "file $outfile cannot be found" unless -e "$outfile";
#	    print "aggregating $outfile\n";
	    my $stranded = ($orientation =~ /F|R/ ? "--stranded" : "");
	    `cat \"$outfile\" | aggregate_scorefilter_vals $stranded >> $filename3`;
	}
	`mv $filename3 $mainResultDir/$output.Unique.cov`;
    }	
    print " done!\n";
}
