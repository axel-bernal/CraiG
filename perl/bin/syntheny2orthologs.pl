#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Gene;
use GenUtils;
use Contig;
use Organism;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $contigFile;
my $dws = 15;
my $ups = 15;
my $usage = "usage : $0 annot1 ctg1 [annot2 ctg2 .. annotN ctgN] < sytheny_file\n";

my @genes = ();
my @contigs = ();
my @organisms = ();
my %contigMappings = ();
if($#ARGV < 0) {
    die $usage;
}

while($#ARGV >= 0) {
    my $org = Organism->new();
    my $annotFile = shift;
    my $contigFile = shift;

    push(@organisms, $org);

    open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
    my %contigs = % {$org->loadContigs(\*CONTIG)};
    push(@contigs, \%contigs);
    close(CONTIG);
    
    open(ANNOT, "$annotFile") || die $usage;
    my %genes = % { $org->loadGenes(\*ANNOT)};
    push(@genes, \%genes);
    close(ANNOT);
}

for(my $i = 0; $i <= $#genes; $i++) {
    foreach my $id (keys %{$genes[$i]}) {
	my $fexon = $genes[$i]->{$id}->{'Exons'}->[0];
	$contigMappings{$fexon->[0].".".$fexon->[1]} = $genes[$i]->{$id};
    }
}

my @orthologs = ();
my %orthoMappings = ();
my $line = '';

while(defined($line = <STDIN>)) {
    my @line = split(/\s+/, $line);
    my $ref = $line[1].".".$line[4];
    my $syn = $line[6].".".$line[9];

    if(!defined($contigMappings{$ref})) {
	next;
    }
    if(!defined($contigMappings{$syn})) {
	next;
    }
    if(defined($orthoMappings{$ref})) {
	if(!defined($orthoMappings{$syn})) {
	    push(@{$orthologs[$orthoMappings{$ref}]}, $contigMappings{$syn});
	    $orthoMappings{$syn} = $orthoMappings{$ref};
	}
    }
    elsif(defined($orthoMappings{$syn})) {
	if(!defined($orthoMappings{$ref})) {
	    push(@{$orthologs[$orthoMappings{$syn}]}, $contigMappings{$ref});
	    $orthoMappings{$ref} = $orthoMappings{$syn};
	}
    }
    else {
	push(@orthologs, [$contigMappings{$ref}, $contigMappings{$syn}]);
	$orthoMappings{$ref} = $orthoMappings{$syn} = $#orthologs;
    }
}

for(my $i = 0; $i <= $#orthologs; $i++) {
    my @ids = @{$orthologs[$i]};
    my $sameBounds = 1;
    my $last_seq = "";
    for(my $j = 0; $j <= $#ids; $j++) {
	my $gene = $ids[$j];
#	my $seq = $gene->{'Seq'};
	my $seq = substr($gene->{'Seq'}, 0, $ups).substr($gene->{'Seq'}, length($gene->{'Seq'}) - $dws + 1, $dws);
	if($last_seq ne "" && $last_seq ne $seq) {
	    $sameBounds = 0;
	    last;
	}
	$last_seq = $seq;
    }

    if(!$sameBounds) {
	next;
    }
    for(my $j = 0; $j <= $#ids; $j++) {
	my $gene = $ids[$j];
	print $gene->{'Id'}, "\t";
	print ($j < $#ids ? "\t" : "\n");
    }
}
