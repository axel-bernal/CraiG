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
my $outContigs = "";
my $jointContigs = 1;
my $gcClasses = "";
my $geoMeans = "";
my $cimm = "";
my $usage = "usage : $0 -gc-classes gc_classes -cimm model -means igenic_geo_mean -outC outContigs -ctg contigs -F locsFile > new_locsFile\nReads configuration from model_file and generate new training data: each sequence containing num_genes genes(1 by default), with artificially generated intergenic (emit_sequence executable must be in the path) in between them and at the flanks. The script also generates coding suppression regions in file outContigs.codr\n";

srand;
&GetOptions("F=s" => \$fastaFile,
            "ctg=s" => \$contigFile,
            "outC=s" => \$outContigs,
	    "geomean=s" => \$geoMeans,
	    "cimm=s" => \$cimm,
            "gc=s" => \$gcClasses,
            "jointContigs=i" => \$jointContigs,
	    );

($fastaFile ne "" && $outContigs ne "") || die "fasta file must be specified. ".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

use enum qw(EXON=0 INTRON=1 INTERGENIC=4);
use enum qw(STRAND_FWD=0 STRAND_COMP=1);
my $org = Organism->new();
my $model = GeneModelProps->new();
my %gcClasses = %{$model->setGCClasses($gcClasses)};

print STDERR "Loading sequences..\n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(FASTA, "<$fastaFile") || die "file $fastaFile not found\n";
open(CODR, ">$outContigs.codr");
open(OUTCONTIGS, ">$outContigs");

my %contigs = % {$org->loadContigs(\*CONTIG)};
my %genes = % { $org->loadGenes(\*FASTA, 'TRAIN', 1)};
&GenUtils::trimStopCodon(\%genes);
my %locsByGCClass = % {$org->getGCClasses('TRAIN', \%gcClasses)};
my $k;
my @geoMeans = split(/,/, $geoMeans);
(scalar(@geoMeans) == scalar keys %gcClasses) || die "#of means != gcClasses\n";
for(my $meanInd = 0; $meanInd < scalar(@geoMeans); $meanInd++) {
    $geoMeans[$meanInd] = ($geoMeans[$meanInd]*$jointContigs)/($jointContigs + 1);
}

my $line;
my $c;
my $j;
my $contig;
my %shifted;
my %entries = %{GenUtils::extractFlankingIGenics(\%contigs, 'TRAIN')};
my @classes = sort {$a <=> $b} keys %entries;
# Generating coding region restrictions. At the beginning all regions are allowed to be coding
my %codingRestrictions = ();
foreach my $cId (keys %contigs) {
    $contig = $contigs{$cId};
    $codingRestrictions{$cId} = "1"x$contig->{'Len'};    
}

for($c = 0; $c <= $#classes; $c++) {
    my $class = $classes[$c];
    #sorting
    my @sorted = @{$entries{$class}};

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
    for(my $i = 0; $i <= $#sorted; $i++) {
 #       print STDERR length($sorted[$i]->[3]),"  $lengths[$i]\n";
        if(length($sorted[$i]->[2]) < $lengths[$i]) {
            my $toAdd = $lengths[$i] - length($sorted[$i]->[2]);
#            print STDERR "emit_sequence -cls $c -M $model->{'IntergenicModel'}.4ctx --length=$toAdd";
            my $seq = `emit_sequence --context=$c --length=$toAdd $cimm`;
            chomp($seq);
#            $seq = lc($seq);
            if($sorted[$i]->[1] eq "5P") {
                $shifted{$sorted[$i]->[0]} = length($seq);
                $sorted[$i]->[2] = $seq.$sorted[$i]->[2];
                $codingRestrictions{$sorted[$i]->[0]} = ("0" x length($seq)).$codingRestrictions{$sorted[$i]->[0]};
            }
            else {
                $sorted[$i]->[2] = $sorted[$i]->[2].$seq;
                $codingRestrictions{$sorted[$i]->[0]} = $codingRestrictions{$sorted[$i]->[0]}.("0" x length($seq));
            }
        }  
    }
# printing
    
    my $geneSeq;
    for(my $i = 0; $i < scalar(@{$entries{$class}}); $i += $jointContigs*2) {
      $contig = $contigs{${$entries{$class}}[$i]->[0]};
	print CODR ">$contig->{'Id'}\n";
	print OUTCONTIGS ">$contig->{'Id'}\n";
      my $minContigs = Utils::min(scalar(@{$entries{$class}}), $i + $jointContigs*2);
      
      for(my $j = $i; $j < $minContigs; $j += 2) {
          $contig = $contigs{${$entries{$class}}[$j]->[0]};
          my $cgenes = $contig->{'Features'}->{'TRAIN'};
          my $min = ${$entries{$class}}[$j]->[3];
          my $max = ${$entries{$class}}[$j + 1]->[3];
          $geneSeq = $contig->getSequence(0, $min, $max, 0, 0, STRAND_FWD);
          
          foreach my $gene (keys %$cgenes) {
              my $exons = $cgenes->{$gene}->{'Exons'};
              my $shift = 0;
              
              if(defined($shifted{${$entries{$class}}[$j]->[0]}))  {
                  $shift = $shifted{${$entries{$class}}[$j]->[0]};
              }
              for(my $e = 0; $e < scalar(@$exons); $e++) {
                  $exons->[$e]->[0] = ${$entries{$class}}[$i]->[0];
                  $exons->[$e]->[1] += $shift;
                  $exons->[$e]->[2] += $shift;
              }
              $cgenes->{$gene}->dump(\*STDOUT, \%contigs, 1);
          }
          print OUTCONTIGS ${$entries{$class}}[$j]->[2], $geneSeq, ${$entries{$class}}[$j + 1]->[2];
          print CODR $codingRestrictions{${$entries{$class}}[$j]->[0]};
          (length(${$entries{$class}}[$j]->[2]) + length($geneSeq) + length(${$entries{$class}}[$j + 1]->[2]) == length($codingRestrictions{${$entries{$class}}[$j]->[0]})) || die "${$entries{$class}}[$i]->[0] ".(length(${$entries{$class}}[$j]->[2]) + length($geneSeq) + length(${$entries{$class}}[$j + 1]->[2]))." ".length($codingRestrictions{${$entries{$class}}[$j]->[0]})."\n";
      }
      print OUTCONTIGS "\n";
      print CODR "\n";
  }
}
close(OUTCONTIGS); 
close(CONTIG);
close(FASTA);
close(CODR);
