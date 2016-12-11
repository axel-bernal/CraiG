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
my $min = -1;
my $max = 0;
my $sources = "";
my $model_prefix = "";
my $gene_locs = "";
my $usage = "usage : $0 -input-template template_file [-sources sources_file] [-mins min_num_sources] [-output-model output_model] [-extract-model-params gene_locations] [-maxs max_num_sources]\n\tBuild the output_model.resources, output_model.filters and output_model.features files that need to be used as input for the craigTrain program. It uses template file to generate these files. The template file should always have a .template extension and all files, including the template, the output model and the sources, should be given as a complete path. Files without complete paths will be assumed to be in the models directory.\n\tWhen option -sources is specified multiple sets of files are generated corresponding to a subset of sources listed in sources_file and defined by -minx and -maxs parameters\n";

&GetOptions("-sources=s" => \$sources,
            "-input-template=s" => \$template,
            "-mins=i" => \$min,
	    "-maxs=i" => \$max,
	    "-output-model=s" => \$model_prefix,
	    "-extract-model-params=s" => \$gene_locs
    );

die "$usage" if $template eq "" || ($model_prefix eq "" && $sources eq "");
my @ids = ();
my @paths = ();

if($sources ne "" && $model_prefix eq "") {
    open(SOURCES, "$sources") || die "Can't open $sources. $usage";
    while($_ = <SOURCES>) {
	$_ =~ /(\S+)\s+(\S+)/;
	push(@ids, $1);
	push(@paths, $2);
    }
    
    close(SOURCES);

    $max = $#paths + 1 if !$max;
    $min = $#paths + 1 if $min < 0;
    $model_prefix = basename($sources);

    if($model_prefix =~ /(\S+)\.([^\.]+)/) {
	$model_prefix = $1;
    }
}
else {
    $min = $max = 0;
}


for(my $i = $min; $i <= $max; $i++) {
    my $mfile = $model_prefix.$i;
    $mfile = $model_prefix if $min == $max;

    open(MODELF, ">$mfile.resources") || die "Can't open $mfile.resources";
    open(TEMPLATE, "$template") || die "Can't open template. $usage";

    while($_ = <TEMPLATE>) {
	if($_ =~ /__RESOURCE_FILE__/) {
	    close(MODELF);
	    open(MODELF, ">$mfile.resources") || die "Can't open $mfile.resources";
	}
	elsif($_ =~ /__FILTER_FILE__/) {
	    close(MODELF);
	    open(MODELF, ">$mfile.filters") || die "Can't open $mfile.features";
	}
	elsif($_ =~ /__FEATURE_FILE__/) {
	    close(MODELF);
	    open(MODELF, ">$mfile.features") || die "Can't open $mfile.features";
	}
	elsif($_ =~ /^__METAINFO__\s+/) {

	    for(my $j = 0; $j < $i; $j++) {
		my $line = $_;
		$line =~ s/__METAINFO__\s+//;
		$line =~ s/__MRID__/$ids[$j]/;
		$line =~ s/__MRPATH__/$paths[$j]/;
		print MODELF $line;
	    }
	}
	elsif($_ =~ /__LIST_MRID__(\S+)/) {
	    my $line = $_;
	    my $detail = $1;
	    if($_ =~ /__N__LIST_MRID__\S+/) { 
		$line =~ s/__N/$i /;
	    }

	    my $list_ids = $ids[0];

	    for(my $j = 1; $j < $i; $j++) {
		$list_ids = $list_ids.$detail." ".$ids[$j];
	    }

	    $line =~ s/__LIST_MRID__/$list_ids/;
	    print MODELF $line;
	}
	else {
	    print MODELF $_; 
	}
    }
    close(TEMPLATE);
    close(MODELF);    

    my $tprefix = basename($template);
    if($tprefix =~ /(\S+)\.([^\.]+)/) {
	$tprefix = $1;
    }
    
    open(RESOURCES, $template);
    my $line;
    while(defined($line = <RESOURCES>)) {
	if($line !~ /^\s*\#/ && $line =~ /\s*\~(\S+)\s*/) {
	    my $aux_resfile = $1;
	    my $cmd = dirname($template)."/".$tprefix.".$aux_resfile $mfile.$aux_resfile";
	    $cmd = $aux_resfile !~ /.top$/ ? "ln -s ".$cmd : "cp ".$cmd;
	    `$cmd` if !-e "$mfile.$aux_resfile";
	}
    }

    close(RESOURCES);

    if($gene_locs ne "") {
	my %cdsStLens = ();
	my %numCPs = ();
	my %cps = ();
	my $partial_fsm = dirname($template)."/".$tprefix.".partial.top";
	my $complete_fsm = dirname($template)."/".$tprefix.".complete.top";
	my $has_lintron = 1 if `grep LONG_INTRON $partial_fsm` ne "";

	HistogramUtils::initGenomicRegLengths(\%cdsStLens, \%numCPs, 64, 64,
					      $has_lintron);

	open(GENES, "$gene_locs") || die "Can't open $gene_locs for reading\n";
	my %genes = % {&GenUtils::readLocsFormat(\*GENES)};
	foreach my $id (keys %genes) {
	    my $gene = $genes{$id};
	    HistogramUtils::extractGenomicRegLengths(\%cdsStLens, 
						     $gene->{'Exons'},
						     $gene->{'Strand'},
						     $has_lintron);
	    
	}
	close(GENES);

	foreach my $keyw (keys %cdsStLens) {
	    my $st_lens = $cdsStLens{$keyw};
	    my $num_cps = $numCPs{$keyw};
	    $cps{$keyw} = HistogramUtils::getLengthControlPoints($st_lens, $num_cps);
	}

	HistogramUtils::createLengthModelParamsFile(\%cps, $partial_fsm,
						    "$mfile.partial.top");
	HistogramUtils::createLengthModelParamsFile(\%cps, $complete_fsm,
						    "$mfile.complete.top");

	my %binning_subs = ();
	HistogramUtils::computeBinningSizes(\%binning_subs, \%cps);
	
	foreach my $keyw (keys %binning_subs) {
	    `sed -i \'s/$keyw/$binning_subs{$keyw}/g\' $mfile.features`;
	}
    }
}

sub computeArrayDeltas {
    my($array, $deltas) = @_;
    for(my $i = 1; $i < scalar(@$array); $i++) {
	push(@$deltas, $array->[$i] - $array->[$i-1]);
    }
}

sub maxElemArray {
    my($lv, $limit) = @_;
    my $max = 0;
    my $max_ind = -1;
    for(my $i = 0; $i < $limit; $i++) {
	if($lv->[$i] > $max) {
	    $max = $lv->[$i];
	    $max_ind = $i;
	}
    }
    return [$max_ind, $max];
}
