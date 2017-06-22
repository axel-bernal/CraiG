#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);


my $usage = "usage : $0 < locations\n";

my $gene;

my %genes = % {&GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};

my %gIds = ();
my %mappedIds = ();
my $i = 0;
for my $id (keys %genes) {
    $gene = $genes{$id};
    my $gId = $gene->{'gId'};
    if(defined($gIds{$gId}) && $gene->contigId() ne $gIds{$gId}) {
        print STDERR "$gId has ", $gene->contigId(), " and $gIds{$gId} as contigs\n";
    }

    if(defined($mappedIds{$gId}->{$gene->contigId()})) {
#        print STDERR "new id for $gene->{'gId'} is ", $mappedIds{$gId}->{$gene->contigId()}, "\n";
        $gene->{'Id'} = $mappedIds{$gId}->{$gene->contigId()}."|".$gene->{'tId'};
    }
    else {
        $gene->{'gId'} = $gene->{'gId'}.".$i";
#        print STDERR "mapped $gId in ", $gene->contigId(), " to ", $gene->{'gId'}."|".$gene->{'tId'}, "\n";
        $gene->{'Id'} = $gene->{'gId'}."|".$gene->{'tId'};
        $mappedIds{$gId}->{$gene->contigId()} = $gene->{'gId'};
        $i++;
    }

    $gIds{$gId} = $gene->contigId();
        
    $gene->dump(\*STDOUT, undef, 1);
}

