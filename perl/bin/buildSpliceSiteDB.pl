#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $usage = "usage : $0 -ctg contigFile  < locations\n";

my $wa = 24;
my $wd = 20;
my $contigFile = "";

&GetOptions("-windowAcc=i" => \$wa,
            "-windowDon=i" => \$wd,
            "-ctg=s" => \$contigFile,
    );


my $org = Organism->new();
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};

foreach my $id (keys %genes) {
    my $exons = $genes{$id}->{'Exons'};
    my $strand = $genes{$id}->{'Strand'};
    my $c = $genes{$id}->{'Contig'};

    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
        if($i == 0) {
            next;
        }
	my $is;
	my $ie;
        if($strand == STRAND_FWD) {    
	    $is = $exons->[$i - 1]->[2] + 1;
            $ie = $exons->[$i]->[1] - 1;
        }
        else {
	    $is = $exons->[$i - 1]->[2] - 1;
	    $ie = $exons->[$i]->[1] + 1;
        }
	my $seq = $c->getSequence($wd/2, $is, $ie, $wa/2);
	my $donseq = substr($seq, 0, $wd);
	my $accseq = substr($seq, length($seq) - $wa, $wa);

	print $donseq," \@\@ ",$accseq, " \@\@ $id \#\# $c->{'Id'} \#\#$i\#$strand\n";
    }
}
