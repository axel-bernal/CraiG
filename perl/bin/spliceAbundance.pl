#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile = "";
my $filter = "gene";
my $truncated = 0;
my $pairs = "";
my $usage = "usage : $0 [-splice-pairs <pairs_file>] [-truncated] -ctg contigFile < locations\nreports splice pairs abundance for each gene's intron\n";

&GetOptions("-ctg=s" => \$contigFile,
            "-filter=s" => \$filter,
	    "-truncated" => \$truncated,
	    "-splice-pairs=s" => \$pairs,
    );
($filter =~ /exon/i || $filter =~ /gene/i) || die "$filter must be either gene or exon\n".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $line;
my %canonDonors = ();
my %canonAcceptors = ();

if(length($pairs)) {
    open(PAIRS, "<$pairs") || die "file $pairs not found\n";
    while($line = <PAIRS>) {
	$line =~ /^(\S\S)\-(\S\S)\s(\d+)/;
	
	if(!defined(%{$canonAcceptors{$1}})) {
	    %{$canonAcceptors{$1}} = ();
	}
	if(!defined(%{$canonDonors{$2}})) {
	    %{$canonDonors{$2}} = ();
	}

	$canonAcceptors{$1}->{$2} = $3;
	$canonDonors{$2}->{$1} = $3;
    }
    close(PAIRS);
}
else {
    %{$canonAcceptors{'GT'}} = ();
    %{$canonDonors{'AG'}} = ();
    $canonAcceptors{'GT'}->{'AG'} = 1;
    $canonDonors{'AG'}->{'GT'} = 1;
}

print STDERR "splices ", join(" ", %canonAcceptors), " ". join(" ",%canonDonors), "\n";
my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my $seq = "";
my %genes = % { $org->loadGenes(\*STDIN)};

OUT:foreach my $id (keys %genes) {
    my $gene = $genes{$id};

    if(scalar(@{$gene->{'Exons'}}) == 1 && !$gene->{'NoStart'} && !$gene->{'NoStop'}) {
        $gene->displayHeader(\*STDOUT, \%contigs);
	print "\tNN-NN-0\n";
        next;
    }

    my $last_don = 'NN';    
    my @abundances = ();
    for(my $i = 0; $i < scalar(@{$gene->{'Exons'}}); $i++) {
        my $exon = $gene->{'Exons'}->[$i];
	
        if($i > 0 || ($i == 0 && $gene->{'NoStart'} == 1 && !$truncated)) {
            $seq = $gene->{'Contig'}->getSequence(2, $exon->[1], $exon->[1], 0, 0, $gene->{'Strand'});            
	    my $acceptor = uc(substr($seq,0,2));
	    my $exp_dons = $canonDonors{$acceptor};
#            if($seq !~ /^ag/i) {# && $seq !~ /^ac/i) {
            if(!defined($exp_dons) || !defined($exp_dons->{$last_don})) {
		push(@abundances, "$last_don-$acceptor-0");
            }
	    else {
		push(@abundances, "$last_don-$acceptor-".$exp_dons->{$last_don});
	    }
        }

        if(($i == scalar(@{$gene->{'Exons'}}) - 1 && $gene->{'NoStop'} == 1 && !$truncated) || $i < scalar(@{$gene->{'Exons'}}) - 1) {
            $seq = $gene->{'Contig'}->getSequence(0, $exon->[2], $exon->[2], 2, 0, $gene->{'Strand'});
	    
	    $last_don = uc(substr($seq, 1, 2));
#            if($seq !~ /gt$/i) { # && $seq !~ /gc$/i && $seq !~ /at$/i)  {
	}
    }
    # if flow reaches here means $gene has only cannonical splice sites
#    print "awk $gene->{'Id'} $gene->{'NoStart'}\n";
    $gene->displayHeader(\*STDOUT, \%contigs);
    print "\t", join(",", @abundances),"\n";
    
}

