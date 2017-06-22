#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile = "";
my $truncated = 0;
my $splice_file = "";
my $canon_pairs = "";
my $usage = "usage : $0 [-canon-pairs <pairs_file>] [-truncated] -ctg contigFile < locations\nReports Splice site pairs present the in the input annotations\n";

&GetOptions("-ctg=s" => \$contigFile,
	    "-truncated" => \$truncated,
	    "-canon-pairs=s" => \$canon_pairs,
    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $line;
my %canonDonors = ();
my %canonAcceptors = ();

if(length($canon_pairs)) {
    open(PAIRS, "<$canon_pairs") || die "file $canon_pairs not found\n";
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

#print STDERR "splices ", join(" ", %canonAcceptors), " ". join(" ",%canonDonors), "\n";
my $org = Organism->new();
my %nonCanonAcc;
my %nonCanonDon;

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my $seq = "";
my %pairs = ();
my %genes = % { $org->loadGenes(\*STDIN)};

OUT:foreach my $id (keys %genes) {
    my $gene = $genes{$id};

    if(scalar(@{$gene->{'Exons'}}) == 1 && !$gene->{'NoStart'} && !$gene->{'NoStop'}) {
        next;
    }

    my $last_don = 'NN';    

    for(my $i = 0; $i < scalar(@{$gene->{'Exons'}}); $i++) {
        my $exon = $gene->{'Exons'}->[$i];
	
        if($i > 0 || ($i == 0 && $gene->{'NoStart'} == 1 && !$truncated)) {
            $seq = $gene->{'Contig'}->getSequence(2, $exon->[1], $exon->[1], 0, 0, $gene->{'Strand'});            
	    my $acceptor = uc(substr($seq,0,2));
#	    print STDERR "pair $last_don -> $acceptor $seq exon $i\n";
	    if($last_don ne 'NN') {
		$pairs{$last_don."-".$acceptor}++;
	    }
        }

        if(($i == scalar(@{$gene->{'Exons'}}) - 1 && $gene->{'NoStop'} == 1 && !$truncated) || $i < scalar(@{$gene->{'Exons'}}) - 1) {
            $seq = $gene->{'Contig'}->getSequence(0, $exon->[2], $exon->[2], 2, 0, $gene->{'Strand'});
	    
	    $last_don = uc(substr($seq, 1, 2));
	}
    }
}

foreach my $pair (keys %pairs) {
	print "$pair $pairs{$pair}\n";
}

