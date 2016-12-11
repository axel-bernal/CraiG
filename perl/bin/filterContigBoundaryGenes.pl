#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile = "";
my $usage = "usage : $0 -ctg contigs < locations > filtered_locations\n\t\tremoves from locations, all those genes that touch the contig\n\t\tboundaries, either at the beginning or end of contig\n";

&GetOptions("-ctg=s" => \$contigFile,
    );

length($contigFile) || die $usage;
my $org = Organism->new();
open(CONTIGS, "$contigFile") || die $usage;
my %contigs = % {$org->loadContigs(\*CONTIGS)};
close(CONTIGS);
my %genes = % { $org->loadGenes(\*STDIN)};

foreach my $id (keys %genes) {
    my $gene = $genes{$id};
    my $clen = $gene->{'Contig'}->{'Len'};
    if($gene->lowestCoord() != 1 && $gene->highestCoord() != $clen) {
	$gene->dump(\*STDOUT, undef, 1);
    }
}
