#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $usage = "usage : $0 -ctg contigFile < locations\n";
# merges $i genes in one single contig
&GetOptions("-ctg=s" => \$contigFile,
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};
open(CONTIG, ">$contigFile.merged") || die "Can't open $contigFile.merged for writing".$usage;
my @cKeys = keys %contigs;
my $counter = 0;
my $genesPerContig = 2;
OUT: while(1) {
    my $base = 0;
    my $id;
    for(my $j = 0; $j < $genesPerContig; $j++) {
	if($counter >= scalar(@cKeys)) {
	    last OUT;
	}
	my $cId = $cKeys[$counter++];
	if($j == 0) {
	    $id = $cId;
	    print CONTIG ">$id\n";
	}
	my $contig = $contigs{$cId};
	foreach my $gId (keys %{$contig->{'Features'}->{'TRAIN'}}) {
	    my $gene = $contig->{'Features'}->{'TRAIN'}->{$gId};
	    my $exons = $gene->{'Exons'};
	    for(my $e = 0; $e < scalar(@$exons); $e++) {
		$exons->[$e]->[0] = $id;
		$exons->[$e]->[1] += $base;
		$exons->[$e]->[2] += $base;
	    }
	    $gene->displayHeader(\*STDOUT, \%contigs);
	    print "\n";
	}
	print CONTIG "$contig->{'Seq'}\n";
	$base += $contig->{'Len'};
    }
}
close(CONTIG);
