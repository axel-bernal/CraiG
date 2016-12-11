#!/usr/bin/perl

#### This script gives sample/experiment specific score (sum of long and short unique reads)  per intron.
### If want to sum up all the scores from different samples/experiments and give one single intron score, use axel_version.pl. There, each intron will be in 1 line only and will have 1 score only. Not able to tell from which sample the score is coming from. 


use strict;
use Getopt::Long;

my $verbose = 0;
my $testNum = 0;
my $returnOverlaps = 0;
my $intronsFile;
my $genesFile; 
my $outputFile = "addGenesOut.txt";
my $printDistinctIntrons = 0;

&GetOptions("intronsFile|i=s" => \$intronsFile,  ## introns file
            "genesFile|g=s" => \$genesFile,  ## genes file
            "outputFile|o=s" => \$outputFile,  ## file for printing output
            "testNum|t=i" => \$testNum,
            "verbose|v!" => \$verbose,
            "printDistinctIntrons|d!" => \$printDistinctIntrons,
            "returnOverlaps!" => \$returnOverlaps,
    );

if(!-e "$intronsFile" || !-e "$genesFile"){
    die "Unable to find either intronsFile '$intronsFile' or genesFile '$genesFile'\n";
}

open(OUT,">$outputFile") || die "Unable to open '$outputFile' for writing\n";

open(F, "$genesFile");
my $prevId;
my @exons;
my %annotIntrons;
my %genes;  ##datastructure ... key sequence source id containing array of genes
my $ctGenes = 0;
foreach my $line (<F>){
    next if $line =~ /^GENE_ID/;
    chomp $line;
    $line =~ s/\"//g;  ##this removes all quotes from the line if they exist
    my @tmp = split("\t",$line);
    if($tmp[0] ne $prevId){ ##first exon of next gene ...
	&processExons(\@exons) if $prevId;
	undef @exons;  ##reset for the next gene
	push(@exons,\@tmp);
	$prevId = $tmp[0]; ##new prevId ...
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
    $data{$tmp[2]}->{$tmp[3]}->{$tmp[4]}->{$tmp[0]}->{$tmp[1]}->{tot} += ($tmp[6] + $tmp[7]);
    push(@{$data{$tmp[2]}->{$tmp[3]}->{$tmp[4]}->{$tmp[0]}->{$tmp[1]}->{introns}},\@tmp);
    $ctLines++;
    print STDERR "Processing line $ctLines\n" if ($ctLines % 10000 == 0 && $verbose);
}
close I;

print STDERR "I've identified ",scalar(keys%data), " sequences containing introns\n";

foreach my $seqid (keys%data){
    foreach my $start (sort{$a <=> $b}keys%{$data{$seqid}}){
	foreach my $end (sort{$a <=> $b}keys%{$data{$seqid}->{$start}}){
	    foreach my $sample (sort {$a <=> $b}keys%{$data{$seqid}->{$start}->{$end}}){
		foreach my $exp (sort {$a <=> $b}keys%{$data{$seqid}->{$start}->{$end}->{$sample}}){
		    ##now lets determine if this is in a gene and if so, is it annotated ...
		    my @iolaps = @{&getGeneIdForIntron($seqid,$start,$end,$sample,$exp)};
		    if($printDistinctIntrons){
			print OUT "$seqid\t$start\t$end\t$sample\t$exp\t$data{$seqid}->{$start}->{$end}->{$sample}->{$exp}->{tot}";
			foreach my $entry (@iolaps) {
			    print OUT "\t$entry->[0]\t$entry->[1]";
			}
			print OUT "\n";
		    }else{
			foreach my $intron (@{$data{$seqid}->{$start}->{$end}->{$sample}->{$exp}->{introns}}){
			    print OUT join("\t",@{$intron});
			    
			    foreach my $entry (@iolaps) {
				print OUT "\t$entry->[0]\t$entry->[1]";
			    }
			    print OUT "\n";
			}
		    }
		}
	    }
	}
    }
}
sub getGeneIdForIntron {
    my($seqid,$start,$end,$sample,$exp) = @_;
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
		push(@iolaps, [$gene->[0], 0]);
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
			push(@iolaps, [$gene->[0],0]);	
		    }
		    else {
			push(@iolaps, [$gene->[0]."(partial_spanning_right)",0]);
		    }
		}
		else {
		    push(@iolaps, [$gene->[0]."(partial_spanning_left)",0]);
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
    my $geneid = $exons[0]->[0];
    my $seqId = $exons[0]->[2];
    my $geneStart =  $exons[0]->[3];
    my $geneEnd =  $exons[-1]->[4];
    push(@{$genes{$seqId}},[$geneid,$geneStart,$geneEnd]);
    for(my $a = 1; $a < scalar(@exons);$a++){  ##start with second exon to generate first intron
	my $start = $exons[$a-1]->[4] + 1;
	my $end = $exons[$a]->[3] - 1;
	$annotIntrons{$seqId}->{$start}->{$end} = $geneid;  ##make id the value for this exon ... all are in genes so non-null
    }
}
