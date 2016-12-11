#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
#### This script sums up all the scores from different samples/experiments and give one single intron score: Each intron will be in 1 line only and will have 1 score only. Not able to tell from which sample the score is coming from. 
### If want to know sample specific score will need to use the script axel_version_ale.pl

use strict;
use Getopt::Long;

my $verbose = 0;
my $testNum = 0;
my $returnOverlaps = 0;
my $intronsFile = "";
my $genesFile = ""; 
my $outputFile = "addGenesOut.txt";
my $printDistinctIntrons = 0;
my $format = "tsv";
&GetOptions("intronsFile|i=s" => \$intronsFile,  ## introns file
            "genesFile|g=s" => \$genesFile,  ## genes file
            "outputFile|o=s" => \$outputFile,  ## file for printing output
            "testNum|t=i" => \$testNum,
            "verbose|v!" => \$verbose,
            "printDistinctIntrons|d!" => \$printDistinctIntrons,
            "returnOverlaps!" => \$returnOverlaps,
	    "genesFormat|f=s" => \$format,
    );

if(!-e "$intronsFile" || !-e "$genesFile"){
    die "Unable to find either intronsFile '$intronsFile' or genesFile '$genesFile'\n";
}

open(OUT,">$outputFile") || die "Unable to open '$outputFile' for writing\n";

open(F, "$genesFile");
my $prevId = "";
my @exons;
my %annotIntrons;
my %genes;  ##datastructure ... key sequence source id containing array of genes
my $ctGenes = 0;
foreach my $line (<F>){
    next if $line =~ /^GENE_ID/;
    chomp $line;
    $line =~ s/\"//g;  ##this removes all quotes from the line if they exist
    my @tmp = split("\t",$line);
    my $geneId = $tmp[0];
    if($format eq "gtf") {
	$tmp[-1] =~ /\S+\s+(\S+)\;\s+\S+\s+\S+\;/;
	$geneId = $1;
    }

    if($geneId ne $prevId){ ##first exon of next gene ...
	&processExons(\@exons) if $prevId;
	undef @exons;  ##reset for the next gene
	push(@exons,\@tmp);
	$prevId = $geneId; ##new prevId ...
    }else{
	push(@exons,\@tmp);
    }
}
&processExons(\@exons) if $prevId; ##have to do the last one
close F;

print STDERR "Processed $ctGenes genes on ",scalar(keys%genes), " sequences\n";

##hmmm ... probably want to sort the genes on each sequence so that are certain that they are in order
## will make it more efficient when scanning to find gene containing an intron
my %sortedGenes;  
print STDERR "Sorting genes\n" if $verbose;
foreach my $seqid (keys%genes){
    @{$sortedGenes{$seqid}} = sort({$a->[1] <=> $b->[1]} @{$genes{$seqid}});
}
print STDERR "Finished sorting genes\n" if $verbose;

## now want to open the introns file and add in gene and whether intron is annotated
## what else should we do here?  should sort by location when we output the file
## could either combine introns or keep sample specific
open(I,"$intronsFile") || die "Unable to open file '$intronsFile'\n";

my %data;
my $ctLines = 0;
while(<I>){
    chomp;
    $_ =~ s/\"//g;  ##removes quotes from line
    my @tmp = split(/\t/,$_);
    $data{$tmp[2]}->{$tmp[3]}->{$tmp[4]}->{tot} += ($tmp[6] + $tmp[7]);
    $data{$tmp[2]}->{$tmp[3]}->{$tmp[4]}->{ct}++;
    push(@{$data{$tmp[2]}->{$tmp[3]}->{$tmp[4]}->{introns}},\@tmp);
    $ctLines++;
    print STDERR "Processing line $ctLines\n" if ($ctLines % 10000 == 0 && $verbose);
}
close I;

print STDERR "I've identified ",scalar(keys%data), " sequences containing introns\n";

foreach my $seqid (keys%data){
#  print STDERR "$seqid has ",scalar(keys%{$data{$seqid}})," starts\n" if $verbose;

    foreach my $start (sort{$a <=> $b}keys%{$data{$seqid}}){
	foreach my $end (sort{$a <=> $b}keys%{$data{$seqid}->{$start}}){
      ##now lets determine if this is in a gene and if so, is it annotated ...
	    my @iolaps = @{&getGeneIdForIntron($seqid,$start,$end)};
	    if($printDistinctIntrons){
		print OUT "$seqid\t$start\t$end\t$data{$seqid}->{$start}->{$end}->{tot}\t$data{$seqid}->{$start}->{$end}->{ct}";
		foreach my $entry (@iolaps) {
		    print OUT "\t$entry->[0]\t$entry->[1]";
		}
		print OUT "\n";
	    }else{
		foreach my $intron (@{$data{$seqid}->{$start}->{$end}->{introns}}){
		    print OUT "$intron->[2]\t$intron->[3]\t$intron->[4]\t$intron->[0]\t$intron->[1]\t$intron->[5]\t$intron->[8]";

		    foreach my $entry (@iolaps) {
			print OUT "\t$entry->[0]\t$entry->[1]";
		    }
		    print OUT "\n";
		}
	    }
	}
    }
}

sub getGeneIdForIntron {
    my($seqid,$start,$end) = @_;
    my $geneid = $annotIntrons{$seqid}->{$start}->{$end};
    if($geneid){
	return [[$geneid,1]];
    }
    
   ##we need to check here that there are any genes on this seqid .. 
    if(! $sortedGenes{$seqid} ){
	return [["seq_not_annotated",0]];
    }

    my @genes = @{$sortedGenes{$seqid}};

  ##now check to see if the intron is before the first or after last
    my $firstGene = $genes[0];
    my $lastGene = $genes[-1];
    if($firstGene->[1] > $end || $lastGene->[2] < $start){
	return [["outside_not_intergenic",0]];
    }
#    if(! $returnOverlaps){ 
#	if($firstGene->[1] < $end && $firstGene->[1] > $start){
#	    return("partial_spanning",0);
#	}
#	if($lastGene->[2] > $start && $lastGene->[2] < $end){
#	    return("partial_spanning",0);
#	}
 #   }
    my @iolaps = ();
    foreach my $gene (@genes){
	if($returnOverlaps) {
	    if($start > $gene->[2] || $gene->[1] > $end) {
		if($gene->[1] > $end) {
		    if($#iolaps >= 0) {
			return \@iolaps;
		    }
		    else {
			return [["intergenic",0]];
		    }
		}
	    }
	    else {
		push(@iolaps, [$gene->[-3], 0]);
	    }
	}else {
	    if($start > $gene->[2] || $gene->[1] > $end) {
		if($gene->[1] > $end) {
		    if($#iolaps >= 0) {
			return \@iolaps;
		    }
		    else {
			return [["intergenic",0]];
		    }
		}
	    }
	    else {
		if($gene->[1] < $start) {
		    if($gene->[2] > $end) {
			push(@iolaps, [$gene->[-3],0]);	
		    }
		    else {
			push(@iolaps, [$gene->[-3]."(partial_spanning_right)",0]);
		    }
		}
		else {
		    push(@iolaps, [$gene->[-3]."(partial_spanning_left)",0]);
		}
	    }
	}
    }
    
    if($#iolaps >= 0) {
	return \@iolaps;
    }
    else {
	return [["intergenic",0]]; ##this will be the case if there are no genes annotated on this seqid}
    }
}

sub processExons {
    my($tmpexons) = @_;
    my @exons = sort {$a->[3] <=> $b->[3]} @{$tmpexons};
    $ctGenes++;
    my $strand = $format eq "gtf" ? $exons[0]->[6] : ($exons[0]->[5] == 0 ? "+" : "-");
    my $geneid = $exons[0]->[0].$strand;
    if($format eq "gtf") {
	$exons[0]->[-1] =~ /\S+\s+(\S+)\;\s+/;
	$geneid = $1;
    }
    my $seqId = $format eq "gtf" ? $exons[0]->[0] : $exons[0]->[2];
    my $geneStart =  $exons[0]->[3];
    my $geneEnd =  $exons[-1]->[4];
    push(@{$genes{$seqId}},[$geneid,$geneStart,$geneEnd]);
#    print "$seqId, $geneid, $geneStart, $geneEnd, $exons[0]->[-1]\n";
    for(my $a = 1; $a < scalar(@exons);$a++){  ##start with second exon to generate first intron
	my $start = $exons[$a-1]->[4] + 1;
	my $end = $exons[$a]->[3] - 1;
	$annotIntrons{$seqId}->{$start}->{$end} = $geneid;  ##make id the value for this exon ... all are in genes so non-null
    }
}
