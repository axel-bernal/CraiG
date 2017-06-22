#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $blockSize = 300;
my $maskPerc = 50;
my $usage = "usage : $0 -mp mask_percentage[50] -bs block_size[300] < igenicRegions";

&GetOptions("-mp=i" => \$maskPerc,
            "-bs=i" => \$blockSize,
	    );

my $org = Organism->new();
my $c;
$/ = "\n>";
while(defined($c = Contig->new('File' => \*STDIN))) {
    my $mask = 0;
    my $blockStart = 0;
    my $blockEnd = 0;
    my $thisBlockSize = $blockSize;
    my $s;
    my $i;
    my @array = ();

    if($c->{'Len'} < $blockSize) {
        $thisBlockSize = $c->{'Len'};
    }    

    for($i = 0; $i < $thisBlockSize; $i++) {
        $s = substr($c->{'Seq'}, $i, 1);
        if($s =~ /[a-z]/) {
            $mask++;
        }
    }

    my $maskContent = ($mask*100.0)/$thisBlockSize > $maskPerc ? 1 : 0;

    for($i = 0; $i < $thisBlockSize/2; $i++) {
        push(@array, $maskContent);
    }
    
    for( ; $i < $c->{'Len'} - $thisBlockSize/2; $i++) {
        $s = substr($c->{'Seq'}, $i - $thisBlockSize/2, 1);
        if($s =~ /[a-z]/) {
            $mask--;
        }
        $s = substr($c->{'Seq'}, $i + $thisBlockSize/2, 1);
        if($s =~ /[a-z]/) {
            $mask++;
        }        
        $maskContent = ($mask*100.0)/$thisBlockSize > $maskPerc ? 1 : 0;
        push(@array, $maskContent);
    }

    for( ; $i <= $c->{'Len'}; $i++) {
        push(@array, $maskContent);
    }

    $blockStart = -1;
    $blockEnd = -1;
    
    for($i = 1; $i <= $#array; $i++) {

        my $offset = 1;
        my $id = $c->{'Id'};

        if($id =~ /^([^\.]+)\.(\d+)/) {
            $id = $1;
            $offset = $2;
        }

        if($blockStart == -1 && $array[$i] == 1) {
            $blockStart = $i + 1;
        }
        
        if($array[$i] == 0 && $array[$i - 1] == 1) {
            $blockEnd = $i + 1;
        }

        if($blockEnd != -1 && $blockStart != -1) {
            if(($blockEnd - $blockStart) < $thisBlockSize) {
                $blockStart = -1;
                $blockEnd = -1;
                next;
            }
            
            my $cString = $c->getSequence(0, $blockStart, $blockEnd, 0, 0, STRAND_FWD);
            print ">$id.", $offset + $blockStart - 1, "\n$cString\n";
            $blockStart = -1;
            $blockEnd = -1;
        }
            
    }
}
$/ = "\n";
