#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use File::Temp;
use strict;
use Getopt::Long;

my $output = "";
my $numContigs = 0;
my $mainResultDir = "";

&GetOptions("-output=s" => \$output,
	    "-mainResultDir=s" => \$mainResultDir,
	    "-numContigs=i" => \$numContigs,
    );

die "Must provide required parameters" unless  $output ne "" &&  $numContigs > 0 && $mainResultDir ne "";

my (undef, $filename1) = File::Temp::tempfile();
my (undef, $filename2) = File::Temp::tempfile();
my $buff_size = $numContigs + 1000;

my @outfiles = `find $mainResultDir -iname \"split.xsignal.*\" -print0 | xargs -0 -s $buff_size ls`;
foreach my $outfile (@outfiles) {
    chomp $outfile;
    die "file $outfile cannot be found" unless -e "$outfile";
    print STDERR "aggregating $outfile\n";
    `cat \"$outfile\" | aggregate_signalfilter_vals --stranded >> $filename1`;
}
`mv $filename1 $mainResultDir/$output.xsignal`;

@outfiles = `find $mainResultDir -iname \"split.xsigscores.*\" -print0 | xargs -0 -s $buff_size ls`;
foreach my $outfile (@outfiles) {
    chomp $outfile;
    die "file $outfile cannot be found" unless -e "$outfile";
    print STDERR "aggregating $outfile\n";
    `cat \"$outfile\" | aggregate_scorefilter_vals --num-cols=3 --static-vals --stranded >> $filename2`;
}

`mv $filename2 $mainResultDir/$output.xsigscores`;
