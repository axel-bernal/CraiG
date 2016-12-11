#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use Gene;
use strict;

use enum qw(STRAND_FWD=0 STRAND_COMP=1);
my $contigFile;
my $usage = "usage : $0 -ctg contigFile < locations > removedLocations\n";

&GetOptions(
            "-ctg=s" => \$contigFile,
            );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};
my %removed = ();

for my $c (keys %contigs) {
    my @genes = sort {&Gene::compareTo($a, $b);} (values %{ $contigs{$c}->{'Features'}->{'TRAIN'}});

    my $gene;
    my $j = 0;

    for(; $j < scalar(@genes); $j++) {
        $gene =  $genes[$j];
        if(scalar(@{$gene->{'Exons'}}) == 1) {
            next;
        }

        if($gene->{'NoStart'} && $j != 0) {
            $removed{$gene->{'Id'}} = 1;
        }
        if($gene->{'NoStop'} && $j != scalar(@genes) - 1) {
            $removed{$gene->{'Id'}} = 1;
        }
    }
}

foreach my $id (keys %genes) {
    if($removed{$id})  {
        next;
    }

    $genes{$id}->dump(\*STDOUT, undef,  1);
}
