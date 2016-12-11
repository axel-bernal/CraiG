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
    my $stop = $gene->{'Contig'}->getSequence(0, $gene->{'Exons'}->[-1]->[2], $gene->{'Exons'}->[-1]->[2], 3, 0, $gene->{'Strand'});
    $stop = substr($stop, 1, 3);
    print "$gene->{'Exons'}->[-1]->[2] $id $stop\n";
    if($stop =~ /taa$/i || $stop =~ /tga$/i || $stop =~ /tag$/i) {
	$gene->{'NoStop'} = 0;

	if($gene->{'Strand'} == STRAND_COMP) {
	    if(!$gene->{'NoStop'}) {
		$gene->{'Exons'}->[-1]->[2] -= 3;
	    }
	}
	else {
	    if(!$gene->{'NoStop'}) {
		$gene->{'Exons'}->[-1]->[2] += 3;
	    }
	}
    }
}

&GenUtils::displayGeneList(\%genes, \*STDOUT, undef, 1);

