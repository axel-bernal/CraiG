#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use GenUtils;
use Gene;
use strict;
use Getopt::Long;
$| = 1;

# -*- perl -*-

my $usage = "usage: $0 [--complete] [--cds-tag=CDS_TAG] [--transcript-tag=TRANSCRIPT_TAG] < genbank file > fasta file\nUse --complete to extract only non-truncated cDNA entries\n";
my $complete = 0;
my $cds_tag = "CDS";
my $transcript_tag = "mRNA";
# loading the information

&GetOptions("--complete" => \$complete,
	    "--cds-tag" => \$cds_tag,
	    "--transcript-tag" => \$transcript_tag,
	    );

my %genes = %{&GenUtils::readGenBankFormat($complete, \*STDIN, 'TEST', $cds_tag, $transcript_tag)};
&GenUtils::writeGTFFormat(\%genes, \*STDOUT);
