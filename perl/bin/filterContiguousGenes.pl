#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use Gene;
use strict;

use enum qw(STRAND_FWD=0 STRAND_COMP=1);
my $contigFile;
my $usage = "usage : $0 -ctg contigFile < locations > newLocations\nFilters genes which are contiguosly occuring over DNA - without intergenic material in the middle\n";

&GetOptions(
            "-ctg=s" => \$contigFile,
            );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*STDIN)};

for my $c (keys %contigs) {
    my @genes = sort {&Gene::compareTo($a, $b);} (values %{ $contigs{$c}->{'Features'}->{'TRAIN'}});
    my $j = 0;
    my $gene;
    for(; $j < scalar(@genes); $j++) {
        $gene =  $genes[$j];
        my $i = $j + 1;
        my $remove = 0;
        while($i < scalar(@genes) && 
              $genes[$i]->lowestCoord() <= $gene->highestCoord() + 1) {
            if($genes[$i]->lowestCoord() == $gene->highestCoord() + 1) {
                print STDERR "$gene->{'Id'} $genes[$i]->{'Id'} transcripts\n";
                if(scalar(@{$gene->{'3UTRExons'}})) {
                    undef(@{$gene->{'3UTRExons'}});
                }
                elsif(scalar(@{$genes[$i]->{'5UTRExons'}})) {
                    undef(@{$genes[$i]->{'5UTRExons'}});
                }
                else {
                    $remove = 1;
                    print STDERR "$gene->{'Id'} removed ", $gene->contigId(), "\n";
                }                
            }
            $i++;
        }
        if(!$remove) {
            $gene->dump(\*STDOUT, undef, 1);
        }
    }
}
