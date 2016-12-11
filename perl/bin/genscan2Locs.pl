#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use GenUtils;
use Contig;
use Gene;
use strict;
use Getopt::Long;
$| = 1;

# -*- perl -*-

my $contigFile = "";
my $preserve_phases = 0;

my $usage = "usage: $0 [--preserve-phases] -ctg contig < locationFile > exonsFile\n";

&GetOptions(
            "--ctg=s" => \$contigFile,
	    "--preserve-phases" => \$preserve_phases,
    );
# loading the information
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %locs = %{&GenUtils::readGenScanFormat(\*STDIN, 'TEST', \%contigs, $preserve_phases)};
&GenUtils::displayGeneList(\%locs, \*STDOUT);
