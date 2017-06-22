#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $ctg = "";
my $outC = "";
my $subsFile = "";
my $contigFile = "";
my $org = Organism->new();

my $usage = "$0 -ctg contigFile -sub subsFile -outC outcontigFile < genes > new_genes\n";
&GetOptions("-ctg=s" => \$contigFile,
            "-outC=s" => \$outC,
            "-sub=s" => \$subsFile,
            );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(SUBS, "<$subsFile") || die "file $subsFile not found\n";
open(OUT_CONTIG, ">$outC") || die "file $outC not found\n";

my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};

my %newContigs = ();
my $line;

while(defined($line = <SUBS>) && $line =~ />?(\S+)\s+(\S+)/) {
    $newContigs{$1} = $2;
}

foreach my $gId (keys %genes) {
    my $gene = $genes{$gId};
    my $contig = $contigs{$gene->contigId()};
    $contig->{'Id'} = $newContigs{$gene->{'Id'}};
    $gene->shiftCoords(0, $contig->{'Len'}, "$contig->{'Id'}");
    $gene->dump(\*STDOUT, undef, 1);
}

foreach my $cKey (keys %contigs) {
    print OUT_CONTIG ">$contigs{$cKey}->{'Id'}\n$contigs{$cKey}->{'Seq'}\n";
}

close OUT_CONTIG;
