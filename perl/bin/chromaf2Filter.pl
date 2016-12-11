#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{HOME}/Applications/Contrib/bioperl-1.5.2_100/blib/lib";
#use lib "$ENV{HOME}/Applications/Contrib/bioperl-1.5.2_100/blib/lib";
use lib "$ENV{CRAIG_HOME}/perl/lib";

use Bio::AlignIO;
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use File::Basename;
use enum qw(STRAND_FWD STRAND_COMP);

my $usage = "$0 -chr chromo [-l] -ctg sub_contigs < maf_alignments\n";

my $chromo = "";
my $subcontig = "";
my $subContigsAreLocs = 0;

&GetOptions("-l" => \$subContigsAreLocs,
            "-ctg=s" => \$subcontig,
	    "-chr=s" => \$chromo,
            );

open(SUBCONTIG, "$subcontig") ||  die $usage;
my $prefix = basename($subcontig);
$prefix =~ /(\S+)\.[^\.]+/;
$prefix = $1;

my $scLengths;
my $contigRanges;
my %seen = ();

if($subContigsAreLocs) {
    my $line;
    ($contigRanges, $scLengths) = @{&GenUtils::loadSubContigLocs(\*SUBCONTIG)};
}
else {
    ($contigRanges, $scLengths) = @{&GenUtils::loadSubContigSeqs(\*SUBCONTIG)};    
}

close(SUBCONTIG);
my $alignio = Bio::AlignIO->new(-fh => \*STDIN, -format => 'maf');
my $last_cid = "";
my $base = 0;
my $range = $contigRanges->{$chromo}->[$base];


while(my $aln = $alignio->seek_aln($range->[2])) {
    my @each_seq = $aln->each_seq();
    my $target = $each_seq[0];
    my ($org, $cid) = split (/\./, $target->id());
    
    ($chromo eq $cid && $target->strand() eq '+' && $aln->is_flush) || die "bad align ".$aln->id()."\n";
    
    if(!defined($contigRanges->{$cid})) {
	$aln = $alignio->next_aln();
        next;
    }

    my $index = $base;
    if($last_cid eq "" || $last_cid ne $cid) {
        $index = $base = 0;
    }
    
#    print $target->id(), " [", $target->start(), " ", $target->end(), "]\t", substr($target->seq(), 0, 20), "\n";

    while($index < scalar(@{$contigRanges->{$cid}})) {
        $range = $contigRanges->{$cid}->[$index];
#        print "\t", $range->[2], " ", $range->[3], "\n";
        if($range->[2] > $target->end()) {
            last;
        }

        if($range->[3] < $target->start()) {
            $index++;
            next;
        }

        my $length =  $range->[3] - $range->[2] + 1;
        my $sub_cid = $range->[0];        
        my @targetSeq = split //, $target->seq();
        
        my $start = Utils::max(1, $target->start() - $range->[2] + 1);
        my $end = Utils::min($length, 
                             $target->end() - $range->[2] + 1);

#        print $target->id(), " - ", $range->[0], " ($start $end)\n";
        my $startI = $range->[2] + $start - $target->start();
        my $endI = $range->[2] + $end - $target->start();
#        print "\tinformant ", $startI, " ", $endI, "\n";
        
        for(my $i = 1; $i < scalar(@each_seq); $i++) {
            my $infSeq = $each_seq[$i];
            my ($inf, undef) = split (/\./, $infSeq->id());
            `touch "$inf.$prefix.cons"`;
            open(FILTER, ">>$inf.$prefix.cons") || die "Can't open $inf.$prefix.cons\n";

            if(!defined($seen{$inf.$sub_cid})) {
                if(! -z "$inf.$prefix.cons") {
                    print FILTER "\n";
                }
                print FILTER ">$sub_cid\t$length\t.\n";
                print FILTER $start - 1, "\t";
            }
            else {
                if($start > $seen{$inf.$sub_cid} + 1) {
                    print FILTER "\n", $start - 1, "\t";
                }
            }

            my $count = 1;
            my @residues = split //, $infSeq->seq();

            for(my $j = 0; $j < @residues; $j++) {
                if ($targetSeq[$j] eq '.' or $targetSeq[$j] eq '-') {
                    next;
                }

                if($count >= $startI && $count <= $endI) {
                    print FILTER $residues[$j];
                }
                elsif($count > $endI) {
                    last;
                }
                $count++;               
            }
            
            $seen{$inf.$sub_cid} = $end;            
            close(FILTER);
        }
                
        $index++;
    }

    if($range->[3] < $target->start()) {
	$base++;
	if($index == scalar(@{$contigRanges->{$cid}})) {
	    last;
	}
    }
    $range = $contigRanges->{$cid}->[$base];    
    $last_cid = $cid;
}

print STDERR "$subcontig done\n";
