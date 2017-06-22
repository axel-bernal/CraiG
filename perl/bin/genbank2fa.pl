#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Contig;
use GenUtils;
use Gene;
use strict;
use Getopt::Long;
$| = 1;
my $usage = "usage: $0 [--complete] [--cds-tag=CDS_TAG] [--transcript-tag=TRANSCRIPT_TAG] < genbank file > fasta file\nUse --complete to extract only non-truncated cDNA entries\n";
# -*- perl -*-
my $complete = 0;
my $cds_tag = "CDS";
my $transcript_tag = "mRNA";

&GetOptions("--complete" => \$complete,
	    "--cds-tag" => \$cds_tag,
	    "--transcript-tag" => \$transcript_tag,
	    );

my $usage = "usage: $0 < genbank file > fasta file\n";
my %contigs = ();
GenUtils::readGenBankFormat(0, \*STDIN, 'TEST', $cds_tag, $transcript_tag, \%contigs, 1);
foreach my $cid (keys  %contigs) {
    my $contig = $contigs{$cid};
    print ">", $contig->{'Id'}, "\n";
    print $contig->{'Seq'}, "\n";
}
