#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Gene;
use GeneModelProps;
use Organism;
use Contig;
use Genefinder;
use Getopt::Long;
use strict;


my $fastaFile = "";
my $contigFile = "";
my $igenicsFile = "";
my $jointContigs = 1;
my $gcClasses = "";
my $geoMeans = "";
my $locs = 0;

my $usage = "usage : $0 -gc gc_classes -geomean igenic_geo_mean [-l] -ctg contigs -igenics intergenic_regions -F locsFile > new_locsFile\nReads configuration from model_file and generates new training data: each sequence containing num_genes genes(1 by default), with added intergenic regions (supplied by igenics parameter) in between them and at the flanks. The script also generates coding suppression regions\n";

srand;
&GetOptions("F=s" => \$fastaFile,
            "ctg=s" => \$contigFile,
	    "geomean=s" => \$geoMeans,
            "gc=s" => \$gcClasses,
	    "igenics=s" => \$igenicsFile,
            "jointContigs=i" => \$jointContigs,
            "l" => \$locs,
	    );

($fastaFile ne "" && $igenicsFile ne "") || die "fasta file must be specified. ".$usage;
($contigFile ne "") || die "Contig File must be specified. ".$usage;

use enum qw(EXON=0 INTRON=1 INTERGENIC=4);
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $org = Organism->new();
my $model = GeneModelProps->new();
my %gcClasses = %{$model->setGCClasses($gcClasses)};

print STDERR "Loading sequences..\n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(FASTA, "<$fastaFile") || die "file $fastaFile not found\n";
open(CODR, ">$contigFile.codr");
open(OUTCONTIGS, ">$contigFile.contigs");

my %contigs = % {$org->loadContigs(\*CONTIG)};
my %genes = % { $org->loadGenes(\*FASTA, 'TRAIN', 1)};

close(CONTIG);
close(FASTA);

&GenUtils::trimStopCodon(\%genes);
my %locsByGCClass = % {$org->getGCClasses('TRAIN', \%gcClasses)};
my $k;
my @geoMeans = split(/,/, $geoMeans);

for(my $meanInd = 0; $meanInd < scalar(@geoMeans); $meanInd++) {
    $geoMeans[$meanInd] = ($geoMeans[$meanInd]*$jointContigs)/($jointContigs + 1);
}

my $line;
my $c;
my $j;
my $contig;
my %segments;
my %entries = %{GenUtils::extractFlankingIGenics(\%contigs, 'TRAIN')};
my @classes = sort {$a <=> $b} keys %entries;
# Generating coding region restrictions. At the beginning all regions are allowed to be coding
my %codingRestrictions = ();
my $cId;
my $offset = 0;

foreach $cId (keys %contigs) {
    $contig = $contigs{$cId};
    $codingRestrictions{$cId} = "1"x$contig->{'Len'};    
}


for($c = 0; $c <= $#classes; $c++) {
    my $class = $classes[$c];
    open(IGENICS, "<$igenicsFile.$class") || die "Cannot open $igenicsFile.$class\n";
    #sorting
    my @sorted = @{$entries{$class}};
    my $igc;
    my $igcSeq = "";
    my $igcId = "";
    my $totalLength = 0;

    for(my $i = 0; $i < $#sorted; $i++) {
        for($j = $i; $j <= $#sorted; $j++) {
            if(length($sorted[$i]->[2]) > length($sorted[$j]->[2])) {
                my $tmp = $sorted[$j];
                $sorted[$j] = $sorted[$i];
                $sorted[$i] = $tmp;
            }            
        }
    }
    
    #generating length distribution
    my @accProbs = @{&Utils::getGeoAccumProbs($geoMeans[$c], $geoMeans[$c]*12)};
    my @lengths = ();
    for(my $i = 0; $i <= $#sorted; $i++) {
        my $length = Utils::pickAtRandom(\@accProbs);
        for($j = 0; $j < $i; $j++) {
            if($length <= $lengths[$j]) {
                last;
            }
        }
#        print "\"",join(" ", @lengths[0 .. $j - 1]), "\" ", $length;
#        print " \"", join(" ", @lengths[$j .. $#lengths]), "\"";
        @lengths = (@lengths[0 .. $j - 1], $length, @lengths[$j .. $#lengths]);
#        print " \| \"",join(" ", @lengths), "\"\n";
    }
    
    #generating sequences
    $/ = "\n>";
    $offset = 0;
    my $igenicOffset = 0;
    for(my $i = 0; $i <= $#sorted; $i++) {
 #       print STDERR length($sorted[$i]->[3]),"  $lengths[$i]\n";
        if(length($sorted[$i]->[2]) < $lengths[$i]) {
            my $toAdd = $lengths[$i] - length($sorted[$i]->[2]);
#            print STDERR "\n$sorted[$i]->[0] class $class to add $toAdd $lengths[$i] ", length($sorted[$i]->[2]), "\n";
            my $seq = "";
            my @seqSegs = ();

            while(1) {
                my $remaining = $toAdd - length($seq);
                if(defined($igc->{'Seq'})) {
                    my $seqRest = length($igc->{'Seq'}) - $offset;
                    if($seqRest < 10) {
                        undef $igc->{'Seq'};
                        next;
                    }
                    my $len2add = ($seqRest > $remaining) ? $remaining : $seqRest;
                    push(@seqSegs,
                         [$igcId, $igenicOffset, $igenicOffset + $len2add - 1]);
                    $seq = $seq.substr($igc->{'Seq'}, $offset, $len2add);
                    $offset += $len2add;
                    $igenicOffset += $len2add;
                    $totalLength += $len2add;
#                    print length($igc->{'Seq'})," $offset $len2add\n";

                    if(length($seq) == $toAdd) {
                        last;
                    }
                } 

                $igc = Contig->new('File' => \*IGENICS);
                if(!defined($igc)) {
                    last;
                }
                if($igc->{'Seq'} =~ /N/i) {
                    $offset = length($igc->{'Seq'}) + 1;
                    next;
                }

                $igcId = $igc->{'Id'};
                $offset = 0;
                $igenicOffset = 0;
                if($igc->{'Id'} =~ /^(\S+)\.(\d+)$/) {
                    $igcId = $1;
                    $igenicOffset = $2 - 1;
                }
            }    

            (length($seq) == $toAdd) || warn "\nInsufficient igenic sequences for class $class ", $igc->{'Id'}, " ", length($igc->{'Seq'})," $toAdd\n";
            chomp($seq);

            if($sorted[$i]->[1] eq "5P") {
                $segments{$sorted[$i]->[0]."5P"} = \@seqSegs;
                if(!$locs) {
                    $sorted[$i]->[2] = $seq.$sorted[$i]->[2];
                }
                $codingRestrictions{$sorted[$i]->[0]} = ("0" x length($seq)).$codingRestrictions{$sorted[$i]->[0]};
            }
            else {
                $segments{$sorted[$i]->[0]."3P"} = \@seqSegs;;
                if(!$locs) {
                    $sorted[$i]->[2] = $sorted[$i]->[2].$seq;
                }
                $codingRestrictions{$sorted[$i]->[0]} = $codingRestrictions{$sorted[$i]->[0]}.("0" x length($seq));
            }
        }  
        print "$class $totalLength\n";
    }
    
    close(IGENICS);
    $/ = "\n";
# printing
    
    my $geneSeq;
    my $entry = $entries{$class};
    for(my $i = 0; $i < scalar(@{$entry}); $i += $jointContigs*2) {
        my $minContigs = Utils::min(scalar(@{$entry}), $i + $jointContigs*2);
        
        for(my $j = $i; $j < $minContigs; $j += 2) {
            $contig = $contigs{$entry->[$j]->[0]};
            $cId = $contig->{'Id'};
            $offset = 1;
            my $cId4disp = $cId;
            if($cId =~ /^(\S+)\.(\d+)$/) {
                $cId = $1;
                $offset = $2;
                $cId4disp = $cId.$offset;
                $contigs{$cId4disp} = $contig;
            }

            if($i == $j) {
                print CODR ">$cId4disp\n";
                print OUTCONTIGS ">$cId4disp", ($locs ? "\t" : "\n");
            }
            else {
                if($locs) { print OUTCONTIGS ";"; }
            }

            my $cgenes = $contig->{'Features'}->{'TRAIN'};
            my $min = $entry->[$j]->[3];
            my $max = $entry->[$j + 1]->[3];
            if(!$locs) {
                $geneSeq = $contig->getSequence(0, $min, $max, 0, 0, STRAND_FWD);
            }
            else {
                $geneSeq = [$cId, $offset, $contig->{'Len'} + $offset - 1];
            }

            my $shift = 0;
            my $k = 0;
            if(defined($segments{$entry->[$j]->[0]."5P"})) {
                for( ; $k < scalar(@{$segments{$entry->[$j]->[0]."5P"}}); $k++) {
                    my $seg = $segments{$entry->[$j]->[0]."5P"}->[$k];
                    $shift += $seg->[2] - $seg->[1] + 1;
                }
            }
            
            foreach my $gene (keys %$cgenes) {
                my $exons = $cgenes->{$gene}->{'Exons'};

                for(my $e = 0; $e < scalar(@$exons); $e++) {
                    $exons->[$e]->[0] = $cId4disp;
                    $exons->[$e]->[1] += $shift;
                    $exons->[$e]->[2] += $shift;
                }
                $cgenes->{$gene}->dump(\*STDOUT, \%contigs, 1);
            }

            if($locs) {
                if(defined($segments{$entry->[$j]->[0]."5P"})) {
                    $k = 0;
                    for( ; $k < scalar(@{$segments{$entry->[$j]->[0]."5P"}}); $k++) {
                        my $seg = $segments{$entry->[$j]->[0]."5P"}->[$k];
                        print OUTCONTIGS "$seg->[0]_",$seg->[1] + 1, "_", 
                        $seg->[2] + 1, ";";
                        
                    }
                }
                print OUTCONTIGS $geneSeq->[0], "_", $geneSeq->[1], "_", $geneSeq->[2];

                if(defined($segments{$entry->[$j]->[0]."3P"})) {
                    $k = 0;
                    for( ; $k < scalar(@{$segments{$entry->[$j]->[0]."3P"}}); $k++) {
                        my $seg = $segments{$entry->[$j]->[0]."3P"}->[$k];
                        print OUTCONTIGS ";$seg->[0]_",$seg->[1] + 1, "_",
                        $seg->[2] + 1;
                    }
                }
            }
            else {
                print OUTCONTIGS $entry->[$j]->[2], $geneSeq, $entry->[$j + 1]->[2];
            }
            print CODR $codingRestrictions{$entry->[$j]->[0]};

            if(!$locs) {
                (length($entry->[$j]->[2]) + length($geneSeq) + length($entry->[$j + 1]->[2]) == length($codingRestrictions{$entry->[$j]->[0]})) || die "$entry->[$i]->[0] ".(length($entry->[$j]->[2]) + length($geneSeq) + length($entry->[$j + 1]->[2]))." ".length($codingRestrictions{$entry->[$j]->[0]})."\n";
            }
        }
        print OUTCONTIGS "\n";
        print CODR "\n";
    }
}
close(OUTCONTIGS); 
close(CODR);
