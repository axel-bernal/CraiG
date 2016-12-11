#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use File::Basename;
use strict;
use HistogramUtils;
use GenUtils;


my $template = "";
my $input_model = "";
my $usage = "usage : $0 -output-template template_file -input-model model_name";

&GetOptions("-input-model=s" => \$input_model,
            "-output-template=s" => \$template);

die $usage if $input_model eq "" || $template eq "";
die $usage if !-e "$input_model.resources" || !-e "$input_model.filters" || !-e "$input_model.features";

my $tprefix = basename("$template");
if($tprefix =~ /(\S+)\.([^\.]+)/) {
    $tprefix = $1;
}

my $line;
open(RESOURCES, "$input_model.resources");
open(OUTPUT, ">$template") || die "Can't open $template for writing";
print OUTPUT "__RESOURCE_FILE__\n";
while(defined($line = <RESOURCES>)) {
    print OUTPUT $line;
    if($line !~ /^\s*\#/ && $line =~ /\s*\~(\S+)\s*/) {
	my $aux_resfile = $1;
	my $ln_cmd = "ln -s $input_model.$aux_resfile ".dirname($template)."/".$tprefix.".$aux_resfile";
	`$ln_cmd` if !-e dirname($template)."/".$tprefix.".$aux_resfile";
    }
}

open(FILTERS, "$input_model.filters");
print OUTPUT "__FILTER_FILE__\n";
while(defined($line = <FILTERS>)) {
    print OUTPUT $line;
}

open(FEATURES, "$input_model.features");
print OUTPUT "__FEATURE_FILE__\n";
while(defined($line = <FEATURES>)) {
    if($line =~ /Feature(\s+)(\S+)Bin4\-(\S+)Length(\s+)(\S+)Bin(\S+)\s\S+\s+\S+\s+\S+\s+(\S+)/) {
	my $state = $3;
	my $offset = $7;
	$line = "Feature$1$2Bin4\-$3Length$4$5Bin$6";

	if($state eq "Intron") {
	    $line .= " null INTRON_NUMBINS INTRON_LOGBASE $offset INTRON_BINSIZE\n";
	}
	elsif($state eq "InterExon") {
	    $line .= " null INTERNAL_EXON_NUMBINS INTERNAL_EXON_LOGBASE $offset INTERNAL_EXON_BINSIZE\n";
	}
	elsif($state eq "SingleExon") {
	    $line .= " null SINGLE_EXON_NUMBINS SINGLE_EXON_LOGBASE $offset SINGLE_EXON_BINSIZE\n";
	}
	elsif($state eq "InitExon") {
	    $line .= " null INIT_EXON_NUMBINS INIT_EXON_LOGBASE $offset INIT_EXON_BINSIZE\n";
	}
	elsif($state eq "TermExon") {
	    $line .= " null LAST_EXON_NUMBINS LAST_EXON_LOGBASE $offset LAST_EXON_BINSIZE\n";
	}
    }
    print OUTPUT $line;
}	

close(FEATURES);
close(OUTPUT);
