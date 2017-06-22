#!/usr/bin/perl -w

## NOTE: this script takes in the full introns file with genes, annotated appended to lines

## probably want to add threshold for min tot ct before consider valid intron: there are 47 different possible samples that could contain a particular ISR.
## also add percentage diff between big and pair average
## could also add in an option for the percent of reads that the less abundant of set must acheive before considering

use strict;
use Getopt::Long;
## might want to use FileHandle;

my $verbose = 0;
my $testNum = 0;
my $intronsFile = "";
my $geneRPKMFile = "";
my $outputFile = "";
my $perSample = 0;

&GetOptions("intronsFile|f=s" => \$intronsFile, 
            "verbose|v!" => \$verbose,  ## if true then prints to stderr info ...
	    "geneRPKM|g=s"=> \$geneRPKMFile,
            "testNum|t=i" => \$testNum,  ## testnumber of lines to read
	    "outputFile|o=s" => \$outputFile,
	    "perSample|p!" => \$perSample, # report per sample
            );


print STDERR "We are running with verbose on\n" if $verbose;

if(!length($intronsFile) || !-e $intronsFile) {
    die "Usage: reportASEvents.pl --intronsFile|f <filename> [--verbose (if true the prints error statements)]\n";
}

open(F, "$intronsFile");

my $ctLines = 0;

my %input = ();
#3 first read in.  Note that putting into the hash structure to compress many instances of an intron into one.
while(my $line = <F>){
    next if $line =~ /^#/;
    next if $line =~ /Chrom\tStart/;
    next if $line =~ /^SOURCE_ID/;
    last if ($testNum > 0 && $ctLines >= $testNum);  ##means if testnum > 0 and $ctLines > $testNum
    $ctLines++; 
    print STDERR "Reading line $ctLines\n" if $ctLines % 10000 == 0;
    chomp $line;
    $line =~ s/\"//g;  ##remove those pesky quotes if they exist
    my @tmp = split(/\t/,$line);
    die "Incorrect intron file format at line \"$line\"\nUse file with sample level information $#tmp\n" if scalar(@tmp) < 7;
    ##need to use the file with all introns including sample info with gene and annotated added.
    #print STDERR "$tmp[2]: $tmp[3] -> $tmp[4]\n" if $verbose;
    my $sample = $perSample ? "$tmp[3].$tmp[4]" : "total";
    push(@{$input{$sample}}, \@tmp);
}
close F;

foreach my $sample (keys %input) {
    open(OUTPUT, ">$outputFile.$sample") || die "Can't open $outputFile.$sample for writing";

    my %datastart;
    for(my $i = 0; $i <= $#{$input{$sample}}; $i++) {
	my $tmp = $input{$sample}->[$i];
	my $seqid = $tmp->[0].$tmp->[6];
	$datastart{$seqid}->{$tmp->[1]}->{$tmp->[2]}->{tot} += ($tmp->[5]);
	$datastart{$seqid}->{$tmp->[1]}->{$tmp->[2]}->{max} = $tmp->[5] if !defined($datastart{$seqid}->{$tmp->[1]}->{$tmp->[2]}->{max}) || $tmp->[5] > $datastart{$seqid}->{$tmp->[1]}->{$tmp->[2]}->{max};
	$datastart{$seqid}->{$tmp->[1]}->{$tmp->[2]}->{ct}++;
	push(@{$datastart{$seqid}->{$tmp->[1]}->{$tmp->[2]}->{introns}},$tmp);
    }
    
## now it is organized by sequence_id and distinct intron ... need to be able to search along sequence by locations for introns
## sort by location and put into simmple structure with array of [start, end]
    my %dataByStart;
    foreach my $seqid (keys%datastart){
	foreach my $start (sort{$a <=> $b}keys%{$datastart{$seqid}}){
	    foreach my $end (sort{$a <=> $b}keys%{$datastart{$seqid}->{$start}}){ 
		push(@{$dataByStart{$seqid}},[$start,$end]);
	    }
	}
    }

##write a block that sorts and loops through the sequences analyzing and applying the rules
    
    foreach my $seqid (sort{$a cmp $b} keys%dataByStart){
	$seqid =~ /(\S+)(\S)/;
	my @introns = @{$dataByStart{$seqid}};  ## make copy to easier to manage the for loop ..
	
	for(my $j = 0; $j<scalar(@introns);$j++){
	    my $bigInt = $introns[$j];
	    my $bigIntData = $datastart{$seqid}->{$bigInt->[0]}->{$bigInt->[1]};
	    my $bigCt = $bigIntData->{tot};
	    next if $bigIntData->{max} < 3 || abs($bigInt->[0] - $bigInt->[1]) > 7500;
	    my $bigAnnot = $bigIntData->{introns}->[0]->[-1];
	    my $geneIds = &getGeneOlapReport($seqid,$bigInt, \%datastart);
	    my $printout = ">$seqid:$bigInt->[0]..$bigInt->[1]\t$bigCt\t$bigAnnot\t$geneIds";
	    
#	print "working with big intron:  $seqid:$bigInt->[0]..$bigInt->[1] ($bigCt,$bigAnnot,$geneIds)\n";
	    # get the list introns which have at least one end lie within $bigInt
	    my @inbounds = sort{$a->[0] <=> $b->[0]} @{&inbound_introns($j, \@introns, $datastart{$seqid})};
	    
#	print "there are ", $#inbounds, " inbound introns for big intron\n";
	    
	    # extract the alternative splicing events as trees in the interval graph of height 1
	    # rooted at bigInt that have no connected leaves -- except through the root
	    my $max_ct = 0;
	    my $max_index = -1;
	    for(my $k = 0; $k <= $#inbounds; $k++) {
		my $int = $inbounds[$k];
		my $intData = $datastart{$seqid}->{$int->[0]}->{$int->[1]};
		my $annot = $datastart{$seqid}->{$int->[0]}->{$int->[1]}->{introns}->[0]->[-1];
		my $ct = $intData->{tot};
		if($max_ct < 0 || $ct > $max_ct) {
		    $max_ct = $ct;
		    $max_index = $k;
		}
	    }
	    
	    my $printout2 = "-1\n";
	    if($max_index >= 0) {
		my $int = $inbounds[$max_index];
		my $intData = $datastart{$seqid}->{$int->[0]}->{$int->[1]};
		$geneIds = &getGeneOlapReport($seqid, $int, \%datastart);
		my $annot = $intData->{introns}->[0]->[-1];
		my $annotLevel = &getAnnotLevel($bigAnnot, $annot, \%datastart);	    
		$printout2 = "$seqid:$int->[0]..$int->[1]\t$max_ct\t$annot\t$geneIds\t$annotLevel\n";
	    }
	    print OUTPUT $printout."\t".$printout2;
	}
    }
    close OUTPUT;
}


sub inbound_introns {
    my($j, $introns, $datastart) = @_;
    my @ret;
    my $start = $introns->[$j]->[0];
    my $end = $introns->[$j]->[1];
    #going backwards 
#    print "inbounds for $start..$end\n";

    for(my $i = 0; $i <= $#$introns; $i++) {
	my $int = $introns->[$i];
	if($start > $int->[1] || $end < $int->[0]) {
	    next if $start > $int->[1];
	    last if $end < $int->[0];
	}
	else {
	    next if $i == $j;
	    my $ct = $datastart->{$int->[0]}->{$int->[1]}->{max};
	    next if $ct < 3 || abs($int->[0] - $int->[1]) > 7500;
#	    print "inbound $int->[0]..$int->[1] fwd $i\n";
	    push(@ret, $int);
	}
    }
    return \@ret;
}

## returns string of gene_ids ... normally will just be one but could be two
sub getGeneOlapReport {
    my($seqid,$int, $datastart) = @_;
    my $data = $datastart->{$seqid}->{$int->[0]}->{$int->[1]}->{introns}->[0];
    my %ret = ();
    for(my $i = 7; $i < $#$data; $i += 2)  {
	$ret{$data->[$i]} = 1;
    }
    return join(",", keys%ret);
}

sub getOlapGeneIds {
    my($seqid,$int, $datastart) = @_;
    my $data = $datastart->{$seqid}->{$int->[0]}->{$int->[1]}->{introns}->[0];
    my %ret = ();
    for(my $i = 7; $i < $#$data; $i += 2)  {
	my $id = $data->[$i];
	if($id =~ /([^\(]+)(\([^\)]+\))*/) {
	    $id = $1;
	    if($1 =~ /^(\S+)(\+|\-)$/) {
		$id = $1;
	    }
	}
       
	$ret{$id} = 1;
    }
    
    return \%ret;
}

# return 0 if 0 annot, 1 if one of pair annot and 2 if big or both l/r annot
sub getAnnotLevel{
    my($bAnnot,$smannots) = @_;
    my $at_least_one_smannot = ($smannots ne "" && $smannots =~ /1/);
    my $all_smannot = ($smannots ne "" && $smannots !~ /0/);
    my $str = "";
#    for(my $i = 0; $i < $#$smannots; $i++) {
#	$str .= $smannots->[$i]."_";
#	$all_smannot &&= ($smannots->[$i] == 1);
#	$at_least_one_smannot = 1 if ($smannots->[$i] == 1);
#    }
    
    if($bAnnot && !$at_least_one_smannot){  ##only big annotated
	return "only big annotated";#.$bAnnot.$all_smannot.\"$smannots\"";
    }
    if(!$bAnnot && $all_smannot) {  ##l & r annotated
	return "all small annotated";#.$bAnnot.$all_smannot.\"$smannots\"";
    }
    if($bAnnot || $at_least_one_smannot){  ## at least one annotated
	return "at least 1 annotated";#.$bAnnot.$all_smannot.\"$smannots\"";
    }
    return 0;
}
