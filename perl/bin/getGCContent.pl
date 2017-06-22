#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $usage = "usage : $0  < fasta file. Computes the GC content of each id in the fasta file\n";
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
use enum qw(GC_CLUSTER_BEG=0 GC_CLUSTER_END=1 GC_CLUSTER_CLASS=2 GC_CLUSTER_SEQ=3);

my $org = Organism->new();
my $window = 1000000;
my $GCClasses = "43,51,57";
my %gcClasses = %{&GenUtils::parseGCClasses($GCClasses)};
my %contigs = % {$org->loadContigs(\*STDIN)};

my %contigClusters = %{&SigUtils::getContigClustersBySignal(\%contigs, \%gcClasses, $window, "G|C", 1)};

foreach my $contig (keys %contigClusters) {
    my $strand;
    for($strand = STRAND_FWD ; $strand <= STRAND_COMP ; $strand++) {
	foreach my $cluster (@{$contigClusters{$contig}->[$strand]}) {
	    print "$contig $strand ", $cluster->[GC_CLUSTER_BEG] - $contigs{$contig}->{'Beg'} + 1, " ", $cluster->[GC_CLUSTER_END] - $contigs{$contig}->{'Beg'} + 1, " $cluster->[GC_CLUSTER_CLASS]\n";
	}
    }
}
