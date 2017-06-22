#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use GenUtils;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $usage = "usage : $0 [-l] -ctg transcripts_info < locs\n";
my $transcriptFile = "";
&GetOptions(
            "-ctg=s" => \$transcriptFile,
    );


(defined($transcriptFile) && length($transcriptFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();
my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my $index;
my $gene;
my $id;
my $geneId;
my $cid;
my $len;
my $offset;
my $exon;
my $i;

open(CONTIGS, "<$transcriptFile") || die "file $transcriptFile not found\n";

my $line;

while($line = <CONTIGS>) {
    my @line = split(/\s+/, $line);
    $id = $line[0];
    $gene = $genes{$id};
    defined($gene) || die $line;
    my $contig = $gene->{'Exons'}->[0]->[0];
#    $contig =~ /MAL(\d+)/;
#    my $n = $1;
#    $line[1] =~ /[^\_]+\_(\d+)/;
    if($line[1] !~ /scf_/) {
	$line[1] =~ /^[^\_]+\_chr(\S+)$/;
	my $m = $1;
#    print STDERR "$contig $m $line[1] $contig\n";
#    if($n != $m || 
	if($contig ne $m || 
	   ($gene->{'Strand'} == STRAND_FWD && $line[4] eq 'Minus') ||
	   ($gene->{'Strand'} == STRAND_COMP && $line[4] eq 'Plus')) {
	    next;
	}
    }

    my $lowestC = $gene->lowestCoord();
    $offset = $line[2] - $lowestC;

    my $thisGene = Gene->new('Gene' => $gene);
    $thisGene->shiftCoords($offset, 1000000000000, $line[1]);
    $thisGene->dump(\*STDOUT, undef, 1);
}

close(CONTIGS);
