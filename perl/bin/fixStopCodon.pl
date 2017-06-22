#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use GenUtils;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "usage : $0 -ctg contigFile < locations\n";

&GetOptions(
            "-ctg=s" => \$contigFile,
    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my %genes = % { $org->loadGenes(\*STDIN)};

for my $id (keys %genes) {
    $gene = $genes{$id};

    if(!$gene->{'NoStop'}) {
	next;
    }
    if(($gene->{'Phase'} + length($gene->{'Seq'})) % 3) {
	next;
    }

# if there is a donor there, maybe we shouldn't move it
    my $new_end = $gene->{'Exons'}->[-1]->[2];
    my $donor = $gene->{'Contig'}->getSequence(0, $new_end, $new_end, 2, 0, $gene->{'Strand'});
    $donor = substr($donor, 1, 2);

    if($donor =~ /gt/i || $donor =~ /gc/) {
	next;
    }

    # looking 3 bp downstream
    if($gene->{'Strand'} == STRAND_COMP) {
	$new_end -= 3;
    }
    else {
	$new_end += 3;
    }

    my $stop = $gene->{'Contig'}->getSequence(0, $new_end, $new_end, 3, 0, $gene->{'Strand'});
    $stop = substr($stop, 1, 3);
#    print "$gene->{'Exons'}->[-1]->[2] $new_end $id $stop\n";

    if($stop =~ /taa$/i || $stop =~ /tga$/i || $stop =~ /tag$/i) {
	$gene->{'NoStop'} = 0;
	$gene->{'Exons'}->[-1]->[2] = $new_end;
    }
}

&GenUtils::displayGeneList(\%genes, \*STDOUT, undef, 1);

