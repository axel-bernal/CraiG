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
my $splice_file = "";
my $canon_pairs = "";
my $usage = "usage : $0 [-output-pairs <pairs_file>] [-canon-pairs <pairs_file>] [-truncated] -ctg contigFile -filter [gene|exon] < locations\nBy default, only bad exons flanked by non-cannonical SSs will be filtered unless they are internal exons or gene is specified as the filter, in which case the whole gene is excluded from the output\n";

&GetOptions("-ctg=s" => \$contigFile,
            "-filter=s" => \$filter,
	    "-truncated" => \$truncated,
	    "-canon-pairs=s" => \$canon_pairs,
	    "-output-pairs=s" => \$splice_file,
    );
($filter =~ /exon/i || $filter =~ /gene/i) || die "$filter must be either gene or exon\n".$usage;
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

    if(scalar(@{$gene->{'Exons'}}) == 1 && !$gene->{'NoStart'} &&
       !$gene->{'NoStop'}) {
        $gene->dump(\*STDOUT);
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
	    my $exp_dons = $canonDonors{$acceptor};

#            if($seq !~ /^ag/i) {# && $seq !~ /^ac/i) {
            if(!defined($exp_dons) || !defined($exp_dons->{$last_don})) {
#		print STDERR "not ok\n";
		if(!defined($exp_dons)) {
		    print STDERR "non canonical acceptor $acceptor @ $gene->{'Id'} ". join("_", @{$exon}[0..2]),"\n";
		}

                $nonCanonAcc{$acceptor}++;
                if($filter =~ /exon/i && ($i == 0 || $i == scalar(@{$gene->{'Exons'}}) - 1) && scalar(@{$gene->{'Exons'}}) > 1) {
                    if($i == 0) {
                        shift(@{$gene->{'Exons'}});
                        $gene->{'NoStart'} = 1;
                        $i--;
                        $gene->computePhase();
                    }
                    else {
                        $gene->{'NoStop'} = 1;
                        pop(@{$gene->{'Exons'}});
                    }
                }
                else {
#                    print "acceptor $id $exon->[1] $exon->[2] $seq $i\n";
                    next OUT;
                }
            }
	}

        if(($i == scalar(@{$gene->{'Exons'}}) - 1 && $gene->{'NoStop'} == 1 && !$truncated) || $i < scalar(@{$gene->{'Exons'}}) - 1) {
            $seq = $gene->{'Contig'}->getSequence(0, $exon->[2], $exon->[2], 2, 0, $gene->{'Strand'});
	    
	    $last_don = uc(substr($seq, 1, 2));
#            if($seq !~ /gt$/i) { # && $seq !~ /gc$/i && $seq !~ /at$/i)  {
            if(!defined($canonAcceptors{$last_don})) {
		print STDERR "non canonical donor $last_don @ $gene->{'Id'} ". join("_", @{$exon}[0..2])," exon $i ", scalar(@{$gene->{'Exons'}}), "\n";
		$nonCanonDon{$last_don}++;
		
                if($filter =~ /exon/i && ($i == 0 || $i == scalar(@{$gene->{'Exons'}}) - 1) && scalar(@{$gene->{'Exons'}}) > 1) {
                    if($i == 0) {
                        shift(@{$gene->{'Exons'}});
                        $i--;
#                        print STDERR "before ", $gene->{'Phase'},"\n";
                        $gene->{'NoStart'} = 1;
                        $gene->computePhase();
#                        print STDERR "after ", $gene->{'Phase'},"\n";
                    }
                    else {
                        $gene->{'NoStop'} = 1;
                        pop(@{$gene->{'Exons'}});
                    }
                }
                else {
#                    print "donor $id $exon->[1] $exon->[2] $seq $i\n";
                    next OUT;
                }
	    }
	    else {
#		print STDERR  "ok\n";
	    }
	}
    }
    # if flow reaches here means $gene has only cannonical splice sites
#    print "awk $gene->{'Id'} $gene->{'NoStart'}\n";
    $gene->dump(\*STDOUT);
}

print STDERR "non-canonical acceptors\n";
foreach my $acc (keys %nonCanonAcc) {
    print STDERR "$acc $nonCanonAcc{$acc}\n";
}

print STDERR "non-canonical donors\n";
foreach my $don (keys %nonCanonDon) {
    print STDERR "$don $nonCanonDon{$don}\n";
}

if(length($splice_file)) {
    open(SPLFD, ">$splice_file") || die "Cannot open $splice_file\n";
    foreach my $pair (keys %pairs) {
	print SPLFD "$pair $pairs{$pair}\n";
    }
    
    close(SPLFD)
}
else {
    print STDERR "splice pairs\n";
    foreach my $pair (keys %pairs) {
	print STDERR "$pair $pairs{$pair}\n";
    }
}
