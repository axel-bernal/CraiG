#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Gene;
use Organism;
use Contig;
use Genefinder;
use Getopt::Long;
use strict;

my $fastaFile = "";
my $contigFile = "";
my $outputFile = "TrainContigScores";
my $GCClasses = "0";
my $usage = "usage : $0 -o outpuFile -ctg contigs -F fastaFile\nOutput in outputFile.signals";

&GetOptions("F=s" => \$fastaFile,
	    "ctg=s" => \$contigFile,
	    "o=s" => \$outputFile,
	    );

($fastaFile ne "") || die "score file must be specified. ".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

use enum qw(EXON=0 INTRON=1 INTERGENIC=4);

my $org = Organism->new();
my $GF = Genefinder->new();
print STDERR "Loading sequences..\n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(FASTA, "<$fastaFile") || die "file $fastaFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*FASTA, 'TRAIN', 1)};

#`rm "$outputFile.signals"`;                                            
print STDERR "Done!\nAnalyzing Contig Signals... \n";                       

foreach my $cKey (keys %{$org->{'Contigs'}}) {                 
    my %contigSignals = ();                                                   
    print STDERR "Done!\nAnalyzing contig ", $cKey, "\n";                
    my $contig = $org->{'Contigs'}->{$cKey};                         
    $contig->getSignals(\%contigSignals, $org->{'SignalPatterns'});    
    $GF->setContigSignals(\%contigSignals);                              
    print STDERR "Done!\nProcessing Signal Information...\n";      
    $GF->classifySignalsWGeneSet($org, 'TRAIN');                            
    print STDERR "Done!\n";                                       
#  training the signals. Output will be put in model/signal                 
#  print STDERR "Training Signals...\n";                    
#  $GF->trainContigSignals($org, \@addSignals, 40, STRAND_FWD); 
#  print STDERR "Done!\nStoring Signals...\n"; 
#  $GF->scoreContigSignals($org, $maxNumExons, $verbose, $outputFile); 
    $GF->storeContigSigScores("$outputFile.$cKey", $org->{'Contigs'});
#    `cat "$outputFile.$cKey.signals" >> "$outputFile.signals"`; 
#    `rm "$outputFile.$cKey.signals"`;     
}                               
close(CONTIG);
close(FASTA);
print STDERR "Done!\n";

