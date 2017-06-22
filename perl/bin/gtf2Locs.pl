#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use GenUtils;
use Getopt::Long;
use strict;
$| = 1;
my $usage = "usage: $0 [--preserve-phases] -ctg contigs -output <prefix output file> [-uniq-id] [-check] [-no_alt] < gff file > locations file\n-check attempts to correct bad annotations, such as phantom small introns (<5), incorrect exon phase, inexistent start, stop, donor or acceptor signals in the sequence.\n- no_alt option (off by default) only outputs the transcript with more exons for each gene\n";
my $contigFile = "";
my $check = 0;
my $no_alt = 0;
my $uniq_id = 0;
my $output = "";
my $preserve_phases = 0;
&GetOptions("--check" => \$check,
            "--ctg=s" => \$contigFile,
            "--no_alt" => \$no_alt,
	    "--uniq-id" => \$uniq_id,
	    "--output=s" => \$output,
	    "--preserve-phases" => \$preserve_phases,
            );

# loading the information
my $file = \*STDIN;
my $i;
my $org = Organism->new();
my $contigs = undef;

my $ofh = select STDERR;
$| = 1;
select $ofh;

if(length($contigFile)) {
    open(CONTIGS, "<$contigFile") || die "Can't open $contigFile\n". $usage."\n";
    $contigs = $org->loadContigs(\*CONTIGS);
}
else {
    warn "Proceeding without contigs. Usage for including contigs is: $usage\n";
}
my %genes = % { &GenUtils::readGFFFormat($file, 'TEST', $contigs, $no_alt, 
					 $uniq_id, $preserve_phases)};

foreach my $source (keys %genes) {
    my $file = \*STDOUT;
    if(length($output)) {
	open(OUT, ">$output.$source.locs") || die "can't open $output.$source.locs\n";
	$file = \*OUT;
    }

    GenUtils::checkAnnotationSanity($genes{$source}, $preserve_phases) if $check;
    &GenUtils::displayGeneList($genes{$source}, $file, undef, 1);	
    close(OUT) if length($output);
}

exit 0;
