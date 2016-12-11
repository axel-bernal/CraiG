#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use File::Basename;
use strict;
use HistogramUtils;
use GenUtils;

my $oparams = "";
my $iparams = "";
my $usage = "usage : $0 -input param_template -output output_param_model < gene_locs.\nUses gene model entries in gene_locs to compute model parameters like length, and bin widths and generate a parameter output model to be used to predict gene_locs-like genes\n";
&GetOptions("-output=s" => \$oparams,
            "-input=s" => \$iparams,
    );

die $usage  if $oparams eq "" || $iparams eq "";

my %cdsStLens = ();
my %numCPs = ();
my %cps = ();
my $has_lintron = 1 if `grep LONG_INTRON $iparams` ne "";

HistogramUtils::initGenomicRegLengths(\%cdsStLens, \%numCPs, 64, 64, 
				      $has_lintron);

my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
foreach my $id (keys %genes) {
    my $gene = $genes{$id};
    HistogramUtils::extractGenomicRegLengths(\%cdsStLens, 
					     $gene->{'Exons'},
					     $gene->{'Strand'},
					     $has_lintron);
}

foreach my $keyw (keys %cdsStLens) {
    my $st_lens = $cdsStLens{$keyw};
    my $num_cps = $numCPs{$keyw};
    $cps{$keyw} = HistogramUtils::getLengthControlPoints($st_lens, $num_cps);
}

HistogramUtils::createLengthModelParamsFile(\%cps, $iparams, $oparams);
