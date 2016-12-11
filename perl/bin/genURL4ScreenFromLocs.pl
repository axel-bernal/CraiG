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

&GetOptions(
            "-ctg=s" => \$contigFile,
    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my $line = <STDIN>;
while(1) {
    my $gene = GenUtils::readLocsFromString($line, undef);
    if(!defined($gene)) {
	last;
    }

    my $begin = $gene->lowestCoord() - 5000;
    $begin = 1 if $begin < 1;
    my $end = $gene->highestCoord() + 5000;
    my $cid = $gene->{'Exons'}->[0]->[0];
    my $c = $contigs{$cid};
    $end = $c->{'Len'} if $c->{'Len'} < $end;
    
    print "http://toxodb.org/cgi-bin/gbrowse/toxodb/?name=$cid:$begin..$end;%1EtgonME49_DBP_Hehl-Grigg_rnaSeq_RSRCCoverage%1EtgonME49_DBP_Hehl-Grigg_rnaSeq_RSRCJunctions%1EGene\n";
    $line = <STDIN>;
}

