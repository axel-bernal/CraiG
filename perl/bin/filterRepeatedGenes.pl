#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < locations\n";

my %genes = ();

my $line = <STDIN>;
my $i;
my $hash_gene;
my $gene;
while(1) {
    $gene = GenUtils::readLocsFromString($line, undef);
    if(!defined($gene)) {
	last;
    }
    
    $hash_gene = $gene->contigId();
    my $exons = $gene->{'5UTRExons'};
    for($i = 0; $i < scalar(@{$exons}); $i++) {
	$hash_gene .= "_".$exons->[$i]->[1]."_".$exons->[$i]->[2];
    }
    $exons = $gene->{'Exons'};
    for($i = 0; $i < scalar(@{$exons}); $i++) {
	$hash_gene .= "_".$exons->[$i]->[1]."_".$exons->[$i]->[2];
    }
    $exons = $gene->{'3UTRExons'};
    for($i = 0; $i < scalar(@{$exons}); $i++) {
	$hash_gene .= "_".$exons->[$i]->[1]."_".$exons->[$i]->[2];
    }

#    print "hash $hash_gene\n";
    my $model_gene = $genes{$hash_gene};
    if(defined($model_gene)) {
	if(defined($exons->[0]->[3])) {
	    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
		my $exon = $exons->[$i];
		my $model_exon = $model_gene->{'Exons'}->[$i];
		if(defined($model_exon->[3])) {
		    $model_exon->[3] += $exon->[3];
		}
		else {
		    push(@{$model_exon}, $exon->[3]);
		}
	    }
	}
    }
    else {
	$genes{$hash_gene} = $gene;
    }
    $line = <STDIN>;
}

while( (undef, $gene) = each(%genes) ) {
    $gene->dump(\*STDOUT, undef, 1);
}
