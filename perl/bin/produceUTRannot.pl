#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $utrs = "";
my $sort = 0;
my $usage = "usage :$0 [-sort] -utr utr_limits < locations\nWhen using the -sort option make sure the file does not contain\nalternatively spliced transcript.. this is not supported yet\n";

&GetOptions("-sort" => \$sort,
            "-utr=s" => \$utrs,
    );

(length($utrs) != 0) || die "utr_limits file not found. ".$usage;

my %locs = % {&GenUtils::readLocsFormat(\*STDIN)};

open(UTRS, "$utrs") || die $usage;
my %utrs = % {&GenUtils::readLocsFormat(\*UTRS)};
close(UTRS);

my $sort_locs;
my $sort_utrs;
if($sort) {
    @{$sort_utrs} = sort {&Gene::compareTo($a, $b);} values %utrs;
    @{$sort_locs} = sort {&Gene::compareTo($a, $b);} values %locs;
}
else {
    @{$sort_locs} = values %locs;
    @{$sort_utrs} = values %utrs;
}

#($#$sort_locs == $#$sort_utrs) || die "number of utrs != number of locs\n";
my $j = 0;
for(my $i = 0; $i <= $#$sort_utrs && $j <= $#$sort_locs; ) {
    my $utr = $sort_utrs->[$i];
    my $loc = $sort_locs->[$j];
    my $b = $utr->contigId() cmp $loc->contigId();
    if(!$b) {
	if($utr->lowestCoord() <= $loc->lowestCoord() &&
	   $utr->highestCoord() >= $loc->highestCoord()) {
	    if($utr->{'Strand'} == $loc->{'Strand'}) {
		my $has_utrs = $loc->addUTRAnnot($utr);
		$sort_locs->[$j]->dump(\*STDOUT);
	    }
	    $i++; $j++;
	}
	elsif(!$sort) {
	    warn "UTR $utr->{'Id'} and CDS $loc->{'Id'} do not align!\nDiscarding\n";
	    $i++; $j++;
	}
	else {
	    if($utr->highestCoord() > $loc->highestCoord()) {
		$j++;
	    }
	    elsif($utr->highestCoord() < $loc->highestCoord()) {
		$i++;
	    }
	}
    }
    elsif($b < 0) {
	print STDERR "pair $utr->{'Id'} $loc->{'Id'} \n";
	die "Your gff3 has problems. Try the --sort option." if !$sort;
	$i++;
    }
    else {
	print STDERR "pair $utr->{'Id'} $loc->{'Id'} \n";
	die "Your gff3 has problems. Try the --sort option" if !$sort;
	$j++;
    }
}


