#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 ilocations < locations.\nRemoves from locations, any gene for which intron boundary signals present in ilocations do not align well with the gene's own intron boundary signals\n";

my $ilocs = shift || die $usage;
open(ILOCS, "$ilocs") || die "Cannot open $ilocs".$usage;

my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my %maps = ();
my $id;
my $xid;
my $mid;
my $grs;
my %discordant_genes = ();
my %gene_ranges = ();

foreach $id (keys %genes) {
    $xid = $genes{$id}->contigId().$genes{$id}->{'Strand'};
    $mid = \$maps{$xid};
    $grs = \$gene_ranges{$xid};

    if(!defined($mid)) {
	%$$mid = ();
	@$$grs = ();
    }
    push(@$$grs, [$id, $genes{$id}->lowestCoord(), $genes{$id}->highestCoord()]);

    my $exons = $genes{$id}->{'Exons'};
    for(my $i = 1; $i < scalar(@$exons); $i++) {
	my $don_pos = 'D'.$exons->[$i - 1]->[2];
	my $acc_pos = 'A'.$exons->[$i]->[1];
	$$mid->{$don_pos} = $id;
	$$mid->{$acc_pos} = $id;
    }
}


my %igenes = % { &GenUtils::readLocsFormat(\*ILOCS, 0, 'TEST', undef)};

foreach $id (keys %igenes) {
    $xid = $igenes{$id}->contigId().$igenes{$id}->{'Strand'};
    $mid = $maps{$xid};
    $grs = $gene_ranges{$xid};

    if(!defined($mid) || !defined($grs)) {
	next;
    }

    my $exons = $igenes{$id}->{'Exons'};
    my $target_gene = "";
    for(my $i = 1; $i < scalar(@$exons); $i++) {
	my $don_pos = $exons->[$i - 1]->[2];
	my $acc_pos = $exons->[$i]->[1];
	my $j = 0;

	for( ; $j < scalar(@$grs); $j++) {
#    	print "$j $grs->[$j]->[1] $don_pos $acc_pos $grs->[$j]->[2]\n";
	    if($grs->[$j]->[1] < $don_pos && $don_pos < $grs->[$j]->[2] &&
	       $grs->[$j]->[1] < $acc_pos && $acc_pos < $grs->[$j]->[2]) {
		last;
	    }
	}

	if(!defined($mid->{'D'.$don_pos}) || !defined($mid->{'A'.$acc_pos})) {
	    # search for the gene
	    if($j < scalar(@$grs)) {
		$discordant_genes{$grs->[$j]->[0]} = 1;
	    }
	}
	elsif($target_gene ne "") {
	    if($j < scalar(@$grs) && $target_gene ne $grs->[$j]->[0]) {
		$discordant_genes{$grs->[$j]->[0]} = 1;
		$discordant_genes{$target_gene} = 1;
	    }
	}
	else {
	    if($j < scalar(@$grs)) {
		$target_gene = $grs->[$j]->[0];
	    }
	}
    }
}

foreach $id (keys %genes) {
    if(!defined($discordant_genes{$id})) {
	$genes{$id}->dump(\*STDOUT, undef, 1);
    }
}
