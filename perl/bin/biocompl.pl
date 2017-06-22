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
my $usage = "usage : $0 < locations >  compl_locations\nIf there are no locations at the beginning or end of the contigs, then the program assumes that locations cid_0_0 and cid_length(cid)_length(cid) are inserted in the input\n";

my %glocs = %{&GenUtils::readLocsFormat(\*STDIN, 0, 'STDIN', undef)};
my @glocs = sort {&Gene::compareTo($a, $b);} (values %glocs);

if(!scalar(@glocs)) {
    return;
}

my $j = 0;
my $glocL = $glocs[0];

for(my $i = 1; $i < scalar(@glocs); $i++) {
    my $cidL =  $glocL->contigId();
    my $glocR = $glocs[$i];
 
    if($glocR->contigId() ne $cidL) {
        $glocL = $glocR;
        next;
    }

    my $cgloc_beg = $glocL->highestCoord() + 1;
    my $cgloc_end = $glocR->lowestCoord() - 1;
    
    my $new_cid = $cidL.".".$cgloc_beg;
    print STDOUT ">$new_cid\t$cidL", "_", $cgloc_beg, "_", $cgloc_end, "\n";
    
    $glocL = $glocR;

}
