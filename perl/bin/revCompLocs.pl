#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile = "";
my $ds = 0;
my $ups = 0;
my $usage = "usage : $0 -ctg contigFile < locations > revcompLocations\n";

&GetOptions(
	    "-ctg=s" => \$contigFile,
	    );

my $org = Organism->new();
my %contigs = ();
(defined($contigFile) && open(CONTIG, "<$contigFile")) || die "contig file $contigFile missing\n".$usage;
%contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};
my %newGenes = %{&GenUtils::revCompLocs(\%genes, 'TRAIN')};
foreach my $id (keys %newGenes) {
    $newGenes{$id}->displayHeader(\*STDOUT, \%contigs);
    print "\n";
    $newGenes{$id}->displaySeq(\*STDOUT);
}
