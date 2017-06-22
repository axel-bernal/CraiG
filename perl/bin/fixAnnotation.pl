#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use GenUtils;
use Getopt::Long;
use strict;
$| = 1;
my $usage = "usage: $0 -ctg contigs < locs file > fixed_locations file\n-correct off-by-one exon coordinate errors\n";
my $contigFile = "";
&GetOptions(
            "-ctg=s" => \$contigFile,
            );

# loading the information
my $file = \*STDIN;
my $i;
my $org = Organism->new();
my $contigs = undef;

if(length($contigFile)) {
    open(CONTIGS, "<$contigFile") || die "Can't open $contigFile\n". $usage."\n";
    $contigs = $org->loadContigs(\*CONTIGS);
}

my %genes = % { $org->loadGenes(\*STDIN)};
&GenUtils::fixAnnotation(\%genes);
foreach my $id (keys %genes) {
    $genes{$id}->dump(\*STDOUT);
}
