#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use GenUtils;
use Getopt::Long;

my $usage = "usage : $0 contig_file < locations > trimmedLocations\n";
my $contigFile = undef;
&GetOptions(
    "-ctg=s" => \$contigFile,
    );

use enum qw(STRAND_FWD STRAND_COMP);
open(CONTIGS,"<$contigFile") || die "Can't open contig file $contigFile\n";
my $org = Organism->new();
my %contigs = %{$org->loadContigs(\*CONTIGS)};
my %genes = %{$org->loadGenes(\*STDIN, 1)};
&GenUtils::trimStopCodon(\%genes);
&GenUtils::displayGeneList(\%genes, \*STDOUT, undef, 1);

