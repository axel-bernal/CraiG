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
my $geneRepFile = "";
my $outputFile = "triples.out";
my $dbOutput = "";
&GetOptions("intronsFile|f=s" => \$intronsFile, 
            "verbose|v!" => \$verbose,  ## if true then prints to stderr info ...
            "outputFile|o=s" => \$outputFile,  ## file for printing output
	    "geneReport|g=s" => \$geneRepFile, ## file for printing gene report
            "testNum|t=i" => \$testNum,  ## testnumber of lines to read
	    "dbOutputFile|dbo=s" => \$dbOutput, ## write down all the name of the databases
            );


print STDERR "We are running with verbose on\n" if $verbose;

if(!length($intronsFile) || !-e $intronsFile || !length($dbOutput)) {
    print STDERR "$intronsFile $dbOutput files not found\n";
    die "Usage: reportASEvents.pl --intronsFile|f <filename> [--dbo dbOutputFile] [--o outputFile] [--verbose (if true the prints error statements)]\n";
}

open(F, "$intronsFile");
##also open a file for writing results
open(O, ">$outputFile");

my %dbNames = ();

my $ctLines = 0;
my %datastart;
my %dataend;
my $dbCount = 0;

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
    my $seqid = $tmp[0].$tmp[6];
    $dataend{$seqid}->{$tmp[2]}->{$tmp[1]}->{tot} += ($tmp[5]);
    $dataend{$seqid}->{$tmp[2]}->{$tmp[1]}->{max} = $tmp[5] if !defined($dataend{$seqid}->{$tmp[2]}->{$tmp[1]}->{max}) || $tmp[5] > $dataend{$seqid}->{$tmp[2]}->{$tmp[1]}->{max};
    $dataend{$seqid}->{$tmp[2]}->{$tmp[1]}->{ct}++;
    push(@{$dataend{$seqid}->{$tmp[2]}->{$tmp[1]}->{introns}},\@tmp);
    $datastart{$seqid}->{$tmp[1]}->{$tmp[2]}->{tot} += ($tmp[5]);
    $datastart{$seqid}->{$tmp[1]}->{$tmp[2]}->{max} = $tmp[5] if !defined($datastart{$seqid}->{$tmp[1]}->{$tmp[2]}->{max}) || $tmp[5] > $datastart{$seqid}->{$tmp[1]}->{$tmp[2]}->{max};
    $datastart{$seqid}->{$tmp[1]}->{$tmp[2]}->{ct}++;
    push(@{$datastart{$seqid}->{$tmp[1]}->{$tmp[2]}->{introns}},\@tmp);
}
close F;

print STDERR "I've parsed $ctLines lines from the file contained in ",scalar(keys%datastart)," sequences\n" if $verbose;

## now it is organized by sequence_id and distinct intron ... need to be able to search along sequence by locations for introns
## sort by location and put into simmple structure with array of [start, end]
my %dataByStart;
foreach my $seqid (keys%datastart){
    foreach my $start (sort{$a <=> $b}keys%{$datastart{$seqid}}){
	foreach my $end (sort{$b <=> $a}keys%{$datastart{$seqid}->{$start}}){  ##sort decending so largest intron first
	    push(@{$dataByStart{$seqid}},[$start,$end]);  ##results in simple nonredundant intron structure
	}
    }
}

my %dataByEndH;
my %dataByEnd;
foreach my $seqid (keys%dataend){
    foreach my $end (sort{$a <=> $b}keys%{$dataend{$seqid}}){
	foreach my $start (sort{$b <=> $a}keys%{$dataend{$seqid}->{$end}}){  ##sort decending so largest intron first
	    push(@{$dataByEnd{$seqid}},[$start,$end]);  ##results in simple nonredundant intron structure
	    %{$dataByEndH{$seqid}} = () if !defined($dataByEndH{$seqid});
	    $dataByEndH{$seqid}->{"$start..$end"} = $#{$dataByEnd{$seqid}}; ##to look up for introns easily
	}
    }
}

##write a block that sorts and loops through the sequences analyzing and applying the rules
my $ctGenesWith = 0;
my $ctGenesWithout = 0;
my $ctTriples = 0;
my $printCounter = 1;
my %geneReport = ();

foreach my $seqid (sort{$a cmp $b} keys%dataByStart){
  my @introns = @{$dataByStart{$seqid}};  ## make copy to easier to manage the for loop ..
  my $dbe = $dataByEnd{$seqid};
  my $dbel = $dataByEndH{$seqid};
  ##now write a for loop as is ordered so need only to look forward in loop ...
#  print "ordered introns(by start) for seqid $seqid\n";
#  for(my $J = 0; $J<scalar(@introns);$J++){
#      print "$introns[$J]->[0]..$introns[$J]->[1]\n";
#  }
#  print "ordered introns(by end) for seqid $seqid\n";
#  for(my $j = 0; $j<= $#$dbe;$j++){
#      print "$dbe->[$j]->[0]..$dbe->[$j]->[1]\n";
#  }
  
  for(my $j = 0; $j<scalar(@introns);$j++){
      my $bigInt = $introns[$j];
      my $bigIntData = $datastart{$seqid}->{$bigInt->[0]}->{$bigInt->[1]};
      my $bigCt = $bigIntData->{tot};
      next if $bigIntData->{max} < 3 || abs($bigInt->[0] - $bigInt->[1]) > 7500;
      my $samples = $bigIntData->{introns};
      my $bigAnnot = $samples->[0]->[-1];
      my %theseDbNames = ();
      my $exper_info = "";
      my $dbName = "";
      my $maxSampleCov = 0;
      for(my $i = 0; $i <= $#$samples; $i++) {
	  $theseDbNames{"$samples->[$i]->[3],$samples->[$i]->[4]"} = 1 if length($dbOutput);
	  $exper_info .= $samples->[$i]->[3].",";
	  $exper_info .= $samples->[$i]->[4].",";
	  $exper_info .= $samples->[$i]->[5];
	  $exper_info .= "," if $i < $#$samples;
	  if($dbName eq "" || $samples->[$i]->[5] > $maxSampleCov) {
	      $maxSampleCov = $samples->[$i]->[5];
	      $dbName = "$samples->[$i]->[3].$samples->[$i]->[4]";
	  }
      }
      $seqid =~ /(\S+)(\S)/;
      my $realid = $1;
      my $realstrand = $2;
      my $biglocation;
      if($realstrand eq "+") {
	  $biglocation = "<".$realid."_".($bigInt->[0] - 1)."_".($bigInt->[0] - 1).";".$realid."_".($bigInt->[1] + 1)."_".($bigInt->[1] + 1).">\t".$bigCt.";".$bigCt;
      }
      else {
	  $biglocation = "<".$realid."_".($bigInt->[1] + 1)."_".($bigInt->[1] + 1).";".$realid."_".($bigInt->[0] - 1)."_".($bigInt->[0] - 1).">\t".$bigCt.";".$bigCt;
      }
      
      my $bigIntGeneIds =  &getOlapGeneIds($seqid,$bigInt);
      my $geneIds = &getGeneOlapReport($seqid,$bigInt);
      $biglocation .= "\t$dbName";
      print "working with big intron:  $seqid:$bigInt->[0]..$bigInt->[1] ($bigCt,$bigAnnot,$exper_info,$geneIds)\n";
      # get the list introns which have at least one end lie within $bigInt
      my @inbounds = sort{$a->[0] <=> $b->[0]} @{&inbound_introns($j, \@introns, $dbe, $dbel, $datastart{$seqid})};
								  
      print "there are ", $#inbounds, " inbound introns for big intron\n";

      # construct the interval graph for this list
      my @ival_graph = @{&introns2interval_graph($bigInt, \@inbounds)};

      # extract the alternative splicing events as trees in the interval graph of height 1
      # rooted at bigInt that have no connected leaves -- except through the root
      my @annot_introns = ($bigAnnot);
      for(my $i = 0; $i <= $#inbounds; $i++) {
	  my $si = $inbounds[$i];
	  my $annot = $datastart{$seqid}->{$si->[0]}->{$si->[1]}->{introns}->[0]->[-1];
	  push(@annot_introns, $annot);
      }

      my @raw_events = @{&find_trees(\@ival_graph, \@annot_introns)};
      print "there are ", $#raw_events, " raw AS events\n";

      #check for phase preservation
      my @alts_events = @{&analyze_phase($bigInt, \@raw_events, \@inbounds)};

      # removing 2-node(intron) trees in which the leaf starts before the root bigInt.
      # those introns already inserted when the leaf was taken as the root
      for(my $k = 0; $k <= $#alts_events; $k++) {
	  next if $#{$alts_events[$k]} != 3;
	  if($inbounds[$alts_events[$k]->[1] - 1]->[0] < $bigInt->[0]) {
	      splice(@alts_events, $k, 1);
	      $k--;
	  }
      }
      print "there are ", $#alts_events, " final AS events after removing invalid ones\n";

      # print
      ##print out triple here ....
      ##should add some methods to comapre total counts of pair vs big, also could check stage specificity

      for(my $i = 0; $i <= $#alts_events; $i++) {
	  my $printout = "$seqid:  $bigInt->[0]..$bigInt->[1] ($bigCt;$bigAnnot;$exper_info;$geneIds) == ";
	  my $smalllocs = "";
	  my $smallannots = "";
	  my $smallweights = "";
	  my %sampleCounts = ();

	  for(my $j = 1; $j < $#{$alts_events[$i]} - 1; $j++) {
	      my $smallInt = $inbounds[$alts_events[$i]->[$j] - 1];
	      my $smallIntData = $datastart{$seqid}->{$smallInt->[0]}->{$smallInt->[1]};
	      my $ct = $smallIntData->{tot};
	      my $exper_info = "";
	      
	      $samples = $smallIntData->{introns};
	      for(my $k = 0; $k <= $#$samples; $k++) {
		  $theseDbNames{"$samples->[$k]->[3],$samples->[$k]->[4]"} = 1
		      if length($dbOutput);
		  $sampleCounts{"$samples->[$k]->[3]\t$samples->[$k]->[4]"} = ["", 0, 0]
		      if !defined($sampleCounts{"$samples->[$k]->[3]\t$samples->[$k]->[4]"});
		  $sampleCounts{"$samples->[$k]->[3]\t$samples->[$k]->[4]"}->[0] = "$samples->[$k]->[3].$samples->[$k]->[4]";
		  $sampleCounts{"$samples->[$k]->[3]\t$samples->[$k]->[4]"}->[1]++;
		  $sampleCounts{"$samples->[$k]->[3]\t$samples->[$k]->[4]"}->[2] += $samples->[$k]->[5];
		  $exper_info .= $samples->[$k]->[3].",";
		  $exper_info .= $samples->[$k]->[4].",";
		  $exper_info .= $samples->[$k]->[5];
		  $exper_info .= "," if $k < $#$samples;
	      }
	      
	      my @sorted_counts = sort {(-1)*($a->[1] != $b->[1])*($a->[1] <=> $b->[1]) + (-1)*($a->[1] == $b->[1])*($a->[2] <=> $b->[2])} (values %sampleCounts);

	      $seqid =~ /(\S+)(\S)/;
	      $realid = $1;
	      $realstrand = $2;
	      if($realstrand eq "+") {
		  $smalllocs = "<".$realid."_".($smallInt->[0] - 1) if $j == 1;
		  $smalllocs .= "_".($smallInt->[0] - 1).";".$realid."_".($smallInt->[1] + 1);
		  $smallweights .= "$ct;";
		  $smalllocs .= "_".($smallInt->[1] + 1).">\t".$smallweights."$ct\t".join("\t", @{$sorted_counts[0]}) if $j == $#{$alts_events[$i]} - 2;
	      }
	      else {
		  $smalllocs = ($smallInt->[0] - 1).">" if $j == 1;
		  $smalllocs = ($smallInt->[1] + 1).";".$realid."_".($smallInt->[0] - 1)."_".$smalllocs;
		  $smallweights = ";$ct".$smallweights;
		  $smalllocs = "<".$realid."_".($smallInt->[1] + 1)."_".$smalllocs."\t$ct".$smallweights."\t".join("\t", @{$sorted_counts[0]}) if $j == $#{$alts_events[$i]} - 2;
	      }

	      if(length($geneRepFile)) {
		  my $smallIntGeneIds = &getOlapGeneIds($seqid,$smallInt);
		  my $ratio = log($bigCt/$ct)/log(10.0);
		  foreach my $oid (keys %$smallIntGeneIds) {
		      next if !defined($bigIntGeneIds->{$oid});
		      @{$geneReport{$oid}} = (0, 0, 0, 0) if !defined($geneReport{$oid});
		      $geneReport{$oid}->[0] = $ratio if abs($ratio) > abs($geneReport{$oid}->[0]);
		      $geneReport{$oid}->[1] += $ratio;
		      $geneReport{$oid}->[2] = ($bigCt + $ct) if ($bigCt + $ct) > $geneReport{$oid}->[2];
		      $geneReport{$oid}->[3]++;
		  }
	      }

	      my $geneIds = &getGeneOlapReport($seqid,$smallInt);
	      my $annot = $smallIntData->{introns}->[0]->[-1];

	      $smallannots .= $annot;
	      $printout .= " -> " if $j > 1;
	      $printout .= "$smallInt->[0]..$smallInt->[1] ($ct;$annot;$exper_info;$geneIds) ";
	  }
	  my $annotLevel = &getAnnotLevel($bigAnnot, $smallannots);

	  if($annotLevel ne 0) {
#        next unless $annotLevel;  ##if uncomment only prints if one of introns is annotated ... should make command line option so don't have to comment in and out
		  ## add tests here for minimum number of reads to consider, how similar the small introns need to be, comparing levels of members of the triple, and ultimately sample specificity.  Make the variables in these command line options so can easily vary and rerun ...
	      if(length($dbOutput)) {
		  foreach my $key (keys %theseDbNames) {
		      $dbNames{$key} = ++$dbCount if !defined($dbNames{$key});
		  }
		  foreach my $key (keys %theseDbNames) {
		      $printout =~ s/$key/$dbNames{$key}/g;
		  }
	      }

	      print O ">$printCounter\t".$printout."\t$alts_events[$i]->[-2]\t$alts_events[$i]->[-1]\t$annotLevel\t$biglocation\t$smalllocs\n";
	      $printCounter++;
	  }
      }
  }
}

close O;

if(length($geneRepFile)) {
    open(GR, ">$geneRepFile") || die "Cannot open $geneRepFile\n";
    foreach my $gid (keys %geneReport) {
	my $idrep = $geneReport{$gid};
	print GR "$gid\t$idrep->[0]\t".
	    ($idrep->[3] ? $idrep->[1]/$idrep->[3] : 0).
	    "\t$idrep->[2]\t$idrep->[3]\n";
	
    }
    close(GR)
}

if(length($dbOutput)) {
    open(DB, ">$dbOutput") || die "Cannot open $dbOutput\n";
    foreach my $key (keys %dbNames) {
	print DB $dbNames{$key}, "\t", $key, "\n";
    }
    close DB;
}

sub analyze_phase {
    my($bi, $asevents, $introns) = @_;
    my @S = sort {(-1)*($#$a <=> $#$b)} @$asevents;
    my @output = ();
    my %seen = ();

    for(my $k = 0; $k <= $#S; $k++) {
#	print "analyzing phase on [", join(",", @{$S[$k]}), "]\n";
	my $all_seen = 1;
        # event_type 0 -> cassette, 1 -> alt 5' 2 -> alt 3' 3 -> alt 3' and 5'       
	my $event_type;
	my $exon_length = 0;
	my $curr_offset = $bi->[0];
	for(my $l = 1; $l <= $#{$S[$k]}; $l++) {	           
	    my $skl = $S[$k]->[$l];
	    my $si = $introns->[$skl - 1];
	    $event_type = ($si->[0] != $bi->[0] ? "alt5" : "cassette") if $l == 1;
	    $exon_length += abs($si->[0] - $curr_offset) - 1;
	    $curr_offset = $si->[1];
	    $all_seen = 0 if !$seen{$skl};
#	    print "cds $si->[0] $si->[1] $bi->[0] $exon_length $all_seen\n";
	}

	my $skl = $S[$k]->[-1];
	if($introns->[$skl - 1]->[1] != $bi->[1]) {
	    if($event_type eq "cassette") {
		$event_type = "alt3";
	    }
	    else {
		$event_type = "alt53";
	    }
	}
	$exon_length += abs($introns->[$skl - 1]->[1] - $bi->[1]);
#	print "final cds $exon_length $all_seen\n";
	my @tmp = @{$S[$k]};
	push(@output, \@tmp);
	push(@{$output[-1]}, $event_type);
	if($exon_length % 3 == 0) {
	    if(!$all_seen) {
		push(@{$output[-1]}, "mod3");
	    }
	    else {
		push(@{$output[-1]}, "mod3seen");
	    }
	    for(my $l = 1; $l <= $#{$S[$k]}; $l++) {	           
		my $skl = $S[$k]->[$l];
#		print "seen $introns->[$skl - 1]->[0] - $introns->[$skl - 1]->[1]\n";
		$seen{$skl} = 1;
	    }
	}
	else {
	    if(!$all_seen) {
		push(@{$output[-1]}, "\!mod3");
	    }
	    else {
		push(@{$output[-1]}, "\!mod3seen");
	    }
	}
    }

    return \@output;
}

sub find_trees {
    my($graph, $annot_introns) = @_;
    my @S = ();
    my @A = ();
    my @N = ();
    my @nolapping = ();
    print "annotated introns :";
    for(my $i = 0; $i <= $#$annot_introns; $i++) {
	print "\tAND\t", $i if  $annot_introns->[$i];
	push(@nolapping, $i) if $annot_introns->[$i];
    }
    print "\n";
    return \@S if $#nolapping < 0;
    for(my $j = 1; $j <= $#{$graph->[0]}; $j++) {
	if($graph->[0]->[$j]) {
	    my $insert_tree = $annot_introns->[0] || $annot_introns->[$j];
	    my @this_nolapping = ();
	    if(!$insert_tree) { # check whether an annotated intron can be potentially added to this tree
		for(my $i = 0; $i <= $#nolapping; $i++) {
		    if(!$graph->[$j]->[$nolapping[$i]]) {
			push(@this_nolapping, $nolapping[$i]);
		    }
		}
		$insert_tree = 1 if $#this_nolapping >= 0;
	    }

	    if($insert_tree) {
		push(@S, [0, $j]);
		push(@A, $annot_introns->[0] || $annot_introns->[$j]);
		push(@N, \@this_nolapping);
		print "insert tree \# $#S [", join(",", @{$S[-1]}), "]\n";
	    }
	}
	
	for(my $k = 0; $k <= $#S; $k++) {
	    my $j_connects_S_k = 0;
	    my $inserted = 0;
	    for(my $l = 1; $l <= $#{$S[$k]}; $l++) {	           
		my $skl = $S[$k]->[$l];
		if($skl == $j) {
		    $inserted = 1;
		    last;
		}
		if($graph->[$j]->[$skl]) {
		    $j_connects_S_k = 1;
		    last;
		}
	    }
	    if(!$inserted && !$j_connects_S_k) {
#		print "no connection between $j and subset S[$k]\n";
		my $insert_tree = $A[$k] || $annot_introns->[$j];
		my @this_nolapping = ();
		if(!$insert_tree) { #check whether an annotated intron can be potentially added to this tree
		    for(my $i = 0; $i <= $#{$N[$k]}; $i++) {
			if(!$graph->[$j]->[$N[$k]->[$i]]) {
			    push(@this_nolapping, $nolapping[$i]);
			}
		    }
		    $insert_tree = 1 if $#this_nolapping >= 0;
		}

		if($insert_tree) {
		    push(@A, $A[$k] || $annot_introns->[$j]);
		    push(@N, \@this_nolapping);
		    my @arr = @{$S[$k]};
		    push(@arr, $j);
		    push(@S, \@arr);
		    print "insert tree \# $#S [", join(",", @arr), "]\n";
		} 
	    }
	}
    }

    return \@S;
}


sub introns2interval_graph {
    my($int, $inbounds) = @_;
    my @ival_graph = ();
    push(@ival_graph, [0]);

    for(my $i = 0; $i <= $#$inbounds; $i++) {
	push(@{$ival_graph[0]}, 1);
    }

    for(my $i = 0; $i <= $#$inbounds; $i++) {
	push(@ival_graph, [1]);
	for(my $j = 0; $j <= $#$inbounds; $j++) {
	    push(@{$ival_graph[$i+1]}, 0);
	}
    }

    for(my $i = 0; $i <= $#$inbounds; $i++) {
	push(@ival_graph, [1]);
	for(my $j = 1; $j <= $#$inbounds; $j++) {
	    push(@{$ival_graph[$i + 1]}, 0);	    
	}
	my @iolaps = ();
	for(my $j = 0; $j <= $#$inbounds; $j++) {
	    next if $j == $i;
#	    print "checking overlap between $i [", $inbounds->[$i]->[0], ",", $inbounds->[$i]->[1], "] and $j [", $inbounds->[$j]->[0],",", $inbounds->[$j]->[1], "]\n";
	    if($inbounds->[$i]->[0] > $inbounds->[$j]->[1] ||
	       $inbounds->[$j]->[0] > $inbounds->[$i]->[1]) {
		last if $inbounds->[$j]->[0] > $inbounds->[$i]->[1];
	    }
	    else {
		if($inbounds->[$j]->[0] < $inbounds->[$i]->[0]) {
		    if($inbounds->[$j]->[1] > $inbounds->[$i]->[1]) {
#			print "found complete\n";
			$ival_graph[$i+1]->[$j+1] = 1; # completely inbound
		    }
		    else {
#			print "found partial spanning right\n";
			$ival_graph[$i+1]->[$j+1] = 1; # partial_spanning_right
		    }
		}
		else {
#		    print "found partial spanning left\n";
		    $ival_graph[$i+1]->[$j+1] = 1; # partial_spanning_left
		}
	    }
	}
    }

    # for(my $i = 0; $i <= $#$inbounds + 1; $i++) {
    # 	for(my $j = 0; $j <= $#$inbounds + 1; $j++) {
    # 	    print $ival_graph[$i]->[$j], "\t";
    # 	}
    # 	print "\n";
    # }

    return \@ival_graph;
}

sub inbound_introns {
    my($j, $introns, $dbe, $dbel, $datastart) = @_;
    my @ret;
    my $start = $introns->[$j]->[0];
    my $end = $introns->[$j]->[1];
    print "inbounds for $start..$end\n";
    #going backwards 
    my $i = $dbel->{"$start..$end"} - 1;
    my %marked = ();
    for( ;$i >= 0; $i--){
	if($start <= $dbe->[$i]->[1]) {  ##has same start or length diff is modulo 3
	    my $ct = $datastart->{$dbe->[$i]->[0]}->{$dbe->[$i]->[1]}->{max};
	    next if $ct < 3 || abs($dbe->[$i]->[0] - $dbe->[$i]->[1]) > 7500;
	    print "inbound $dbe->[$i]->[0]..$dbe->[$i]->[1] rev $i\n";
	    push(@ret,$dbe->[$i]);
	    $marked{"$dbe->[$i]->[0]..$dbe->[$i]->[1]"} = 1;
	} else{
	    last;  ##since are sorted if the value is different then we've gone beyond this start site
	}
    }
    # going forward
    $i = $j + 1;  ##start with next one
    for( ;$i <= $#$introns;$i++){
	if($end >= $introns->[$i]->[0]) {  ##has same start or length diff is modulo 3	    
	    if(!$marked{"$introns->[$i]->[0]..$introns->[$i]->[1]"}) {
		my $ct = $datastart->{$introns->[$i]->[0]}->{$introns->[$i]->[1]}->{max};
		next if $ct < 3 || abs($introns->[$i]->[0] - $introns->[$i]->[1]) > 7500;
		print "inbound $introns->[$i]->[0]..$introns->[$i]->[1] fwd $i\n";
		push(@ret, $introns->[$i]);
	    }
	}
	else{
	    last;  ##since are sorted if the value is different then we've gone beyond this start site
	}
    }
    return \@ret;
}

## returns string of gene_ids ... normally will just be one but could be two
sub getGeneOlapReport {
    my($seqid,$int) = @_;
    my $data = $datastart{$seqid}->{$int->[0]}->{$int->[1]}->{introns}->[0];
    my %ret = ();
    for(my $i = 7; $i < $#$data; $i += 2)  {
	$ret{$data->[$i]} = 1;
    }
    return join(",", keys%ret);
}

sub getOlapGeneIds {
    my($seqid,$int) = @_;
    my $data = $datastart{$seqid}->{$int->[0]}->{$int->[1]}->{introns}->[0];
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
