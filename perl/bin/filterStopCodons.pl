#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile = "";
my $usage = "usage : $0 -ctg contigFile - < locations";

&GetOptions("-ctg=s" => \$contigFile,
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my $seq;
my %genes = % { $org->loadGenes(\*STDIN)};

foreach my $id (keys %genes) {
    my $gene = $genes{$id};
    if($gene->computePhase()) {
        $gene->dump(\*STDOUT);
        next;        
    }
}
