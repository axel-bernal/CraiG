#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use strict;

use File::Basename;

my $usage = "$0 -ctg sub_contigs -mafdir maf_alignments_subdir -out out_dir\n";

my $subcontig = "";
my $mafDir = 0;
my $outDir = "";
&GetOptions("-mafDir=s" => \$mafDir,
            "-ctg=s" => \$subcontig,
	    "-out=s" => \$outDir,
            );

my $prefix = basename($subcontig);
$prefix =~ /(\S+)\.[^\.]+/;
$prefix = $1;

open(SUBCONTIG, "$subcontig") ||  die $usage;

#my $qsub = "nohup ";
my $qsub = "qsub -p 0 -j y -o nohup.out -q all.q -cwd -b yes -V ";
my  $prg = "$ENV{CRAIG_HOME}/perl/bin/chromaf2Filter.pl";
my $command;

my $scLengths;
my $contigRanges;

($contigRanges, $scLengths) = @{&GenUtils::loadSubContigLocs(\*SUBCONTIG)};

close(SUBCONTIG);

foreach my $cid (keys %$contigRanges) {
    my $compressed = 1;
    my $maf = $mafDir."/$cid.maf.gz";
    if(!-e $maf) {
	$compressed = 0;
	$maf = $mafDir."/$cid.maf";
    }
    if(!-e $maf) {
	warn "$maf cannot be found\n";
	exit;
    }
    
    if(!-e "$outDir/$prefix") {
	`mkdir $outDir/$prefix`;
    }

    my $index = 0;
    my $output = "";
    my $fileno = 1;
    while($index < scalar(@{$contigRanges->{$cid}})) {
	my $range = $contigRanges->{$cid}->[$index];
	$output .= ">$range->[0]\t$cid"."_"."$range->[2]"."_"."$range->[3]\n";
	
	if($index && !($index % 100)) {
	    my $subdir = "$outDir/$prefix/$cid.$fileno";
	    if(-e "$subdir") {
		`rm -rf $subdir`;
	    }
	    `mkdir $subdir`;
	    my $filename = "$subdir/$prefix.subcontigs";
	    open(TMP, ">$filename");
	    print TMP $output;
	    close TMP;
	    $command = "zcat $maf | ".$prg;

	    if(!$compressed) {
		$command = "cat $maf | ".$prg;
	    }

	    $command .= " -chr $cid -l -ctg $prefix.subcontigs";
	    my $prg_file = "$subdir/$cid.$fileno.sh";
	    open(TMP, ">$prg_file");
	    print TMP $command;
	    close(TMP);
	    `chmod +x $prg_file`;
	    `cd $subdir && $qsub $prg_file`;
	    
	    $output = "";
	    $fileno++;
	}
	$index++;
    }

    if(length($output)) {
	my $subdir = "$outDir/$prefix/$cid.$fileno";
	if(-e "$subdir") {
	    `rm -rf $subdir`;
	}
	`mkdir $subdir`;
	my $filename = "$subdir/$prefix.subcontigs";
	open(TMP, ">$filename");
	print TMP $output;
	close TMP; 

	$command = "zcat $maf | ".$prg;
	
	if(!$compressed) {
	    $command = "cat $maf | ".$prg;
	}
	
	$command .= " -chr $cid -l -ctg $prefix.subcontigs";
	my $prg_file = "$subdir/$cid.$fileno.sh";
	open(TMP, ">$prg_file");
	print TMP $command;
	close(TMP);
	`chmod +x $prg_file`;
	`cd $subdir && $qsub $prg_file`;
	
    }

}

