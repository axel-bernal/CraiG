#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
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
foreach my $id (keys %genes) {
    my $contig = $genes{$id}->{'Contig'};
    for(my $i = 0; $i < scalar(@{$genes{$id}->{'Exons'}}); $i++) {
        my $exon = $genes{$id}->{'Exons'}->[$i];
        my $exonSeq = $contig->getSequence(0, $exon->[1], $exon->[2], 0, 0, $genes{$id}->{'Strand'});
        my $Ns = 0.0;
        for(my $j = 0; $j < length($exonSeq); $j++) {
            if(substr($exonSeq, $j, 1) eq 'N') {
                $Ns++;
            }
        }
#        print $exonSeq," N ", $Ns, "\n";
        if($Ns/length($exonSeq) > 0.5) {
            splice(@{$genes{$id}->{'Exons'}}, $i, 1);
            $i--;
        }
    }
    if(scalar(@{$genes{$id}->{'Exons'}}) > 0) { 
        $genes{$id}->dump(\*STDOUT, undef, 1);
    }
}

