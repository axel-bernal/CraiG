#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use GenUtils;

$| = 1;

# -*- perl -*-

$usage = "usage: $0 < gff file > locations file\n";
# loading the information
my $file = \*STDIN;
my $i;
my %genes = % { &GenUtils::readGFF2Format($file, "cds")};
&GenUtils::displayGeneList(\%genes, \*STDOUT);
for($j = 0; $j < scalar(@genes); $j++) {
    my  $gene = $genes[$j];
    print ">$gene->[0]->[0]_$j ";
    for($i = 0; $i < scalar(@$gene); $i++) {
	$feature = $gene->[$i];
	if(!($i % 2)) { # exon
	    $feature->[1] eq "exon" || die @$feature;
	    print "$feature->[0]_$feature->[2]_$feature->[3]";
	    if($i < scalar(@$gene) - 1) {
		print " ; ";
	    }
	    else {
		print "\n";
	    }
	}
	else {
	    $feature->[1] eq "intron" || die @$feature;	
	}
    }
}
