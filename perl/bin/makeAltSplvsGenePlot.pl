#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use File::Basename;
use enum qw(STRAND_FWD STRAND_COMP);

# for matching cds against an annotation, the annotation must appear 
# as the first member in all clusters if at all
my $annotation = "";
my $subcontigs = "";
my %ratios = ();
my @abund = ({},{},{},{});
my @ratios = ({},{},{},{});
my @abund_ratios = ({},{},{},{});

my $usage = "usage : $0 -locs annotation -sc subcontigs < altspl_report\naltspl_report is produced with tool report_rnaseq_asevents\n";
&GetOptions("-locs=s" => \$annotation,
	    "-sc=s" => \$subcontigs);

$annotation ne "" || die $usage;

my $id;
my $i;
my $gene;
my $olap_gene;

my %subcontigs = ();
if($subcontigs ne "") {
    open(SUBCONTIGS, "<$subcontigs") || die "file $subcontigs not found\n";
    %subcontigs = % { &GenUtils::readLocsFormat(\*SUBCONTIGS, 0, 'CONTIGS', undef)};
    close(SUBCONTIGS);
}

open(LOCS, "<$annotation") || die "Can't open $annotation\n".$usage;
my %genes = % { &GenUtils::readLocsFormat(\*LOCS, 0, 'TRAIN', undef)};
close(LOCS);

if($subcontigs ne "") {
    foreach $id (keys %genes) {
	$genes{$id}->moveLocation2OrigContigs(\%subcontigs);
    }
}

my @genes = sort {&Gene::compareTo($a, $b);} (values %genes);
my %genes = ();
my $line;

while(defined($line = <STDIN>) && $line !~ /^stats\[/) { ; }

while(defined($line) && $line =~ /^stats\[(\d+)\]\[(\d+)\]\[(\d+)\]\s+\=\s+(\d+)\s+(\S.+)$/) {
    my $type = $1;
    my $abundance = $2;
    my $rt = $3;
    my $real_measures = $5;
    my @real_measures = split(/\s+/, $real_measures);
#    print "$line $type $abundance $rt\n";
    my $rgenes = GenUtils::readLocsFormat2Array(\*STDIN, 0, 'OLAP', \$line);
    my $num_events = scalar(@$rgenes);
    if($subcontigs ne "") {
	for($i = 0; $i < $num_events; $i++) {
	    $gene = $rgenes[$i]->[1];
	    $gene->moveLocation2OrigContigs(\%subcontigs);
	}
    }

    $abund[$type]->{$abundance} += $num_events;
    $abund[3]->{$abundance} += $num_events;
    $ratios[$type]->{$rt} += $num_events;
    $ratios[3]->{$rt} += $num_events;
    
    if(!defined($abund_ratios[$type]->{$abundance})) {
	%{$abund_ratios[$type]->{$abundance}} = ();
    }
    $abund_ratios[$type]->{$abundance}->{$rt} += $num_events;
    
    if(!defined($abund_ratios[3]->{$abundance})) {
	%{$abund_ratios[3]->{$abundance}} = ();
    }
    $abund_ratios[3]->{$abundance}->{$rt} += $num_events;
# marginalize abundance
    if(!defined($ratios{$rt})) {
	@{$ratios{$rt}} = ()
    }    
    
#detailed table of events
    print "Per-Event Table\n";
    for($i = 0; $i < $num_events; $i++) {
	my $rabund = $real_measures[2*$i];
	my $rratio = $real_measures[2*$i + 1];	    
	$gene = $rgenes->[$i];
	push(@{$ratios{$rt}}, $gene->[1]);
	my $exons = $gene->[1]->{'Exons'};
	my %olapGenes = ();
#    print $rt, "-> ", scalar(@{$ratios{$rt}}), " ", scalar(@genes), "\n";
	GenUtils::computeGeneOverlaps(\@genes,[$gene->[1]],\%olapGenes,0,0,0);
	
	foreach my $index (keys %olapGenes) {
	    $olap_gene = $genes[$index];
	    my $oid = $olap_gene->{'tId'};
	    if(!defined($genes{$oid})) {
		$genes{$oid} = [$#{$olap_gene->{'Exons'}}, $rabund, $type, $rabund, $rratio, $type, altplsScore($type, $abundance, $rt)];
	    }
	    else {
		if($genes{$oid}->[1] < $rabund) {
		    $genes{$oid}->[1] = $rabund;
		    $genes{$oid}->[2] = $type;
		}
		if($type != 3) {
		    my $score = altSplScore($type, $abundance, $rt);
		    if($genes{$oid}->[6] < $score) {
			$genes{$oid}->[4] = $rabund;
			$genes{$oid}->[5] = $rratio;
			$genes{$oid}->[6] = $score;
		    }
		}
	    }
	}

	($rt >= 4 || $abundance >= 3) || next;

	if(scalar(keys %olapGenes)) {
	    foreach my $index (keys %olapGenes) {
		$olap_gene = $genes[$index];
		print "GeneId\t\t", $olap_gene->{'tId'}, "\n";
		print "Strand\t\t", ($olap_gene->{'Strand'} == STRAND_FWD ? "Sense" : "Antisense"), "\n";
		print "MaxCDSposition\t",$olap_gene->{'Exons'}->[0]->[0],":", $olap_gene->lowestCoord(), "..", $olap_gene->highestCoord(), "\n";
		print "Gene Model\n";
		$olap_gene->displayHeader(\*STDOUT, undef, ",", ":", "..");
		print "\n";
	    }
	}
	else {
	    print "GeneId\t\tN.A. (no overlap to known gene\n";
	    print "Strand\t\t", ($gene->[1]->{'Strand'} == STRAND_FWD ? "Sense" : "Antisense"), "\n";
	    print "MaxCDSposition\t",$exons->[0]->[0],":", $gene->[1]->lowestCoord(), "..", $gene->[1]->highestCoord(), "\n";
	    print "Model\tN.A.\n";
	}

	print "Abundance Main Alternative Form(log)\t", $abundance, "\n";
	print "Ratio\t\t", $rt*5, "\n";
	print "Type\t", ($type == 2 ? "Exon Skip" : ($type == 1 ? "Multiple Overlap of Introns" : "Alternative 3 or 5 splice end")), "\n";
	print "Main Alternative Form\n\t";
	print $exons->[0]->[0],":", $exons->[0]->[1],"..", $exons->[0]->[2], "\n";
	print "Secondary Alternative Form\n\t";
	for(my $j = 1; $j < scalar(@$exons); $j++) {
	    print $exons->[$j]->[0],":", $exons->[$j]->[1],"..", $exons->[$j]->[2];
	    if($j < scalar(@$exons) - 1) {
		print ",";
	    }
	}	
	print "\n";
    }
}

#detailed table of genes
print "Per-Gene Table\n";
my @gene_abunds = ();
foreach my $id (keys %genes) {
    $gene = $genes{$id};
    $gene_abunds = 
    print join("\t", @$gene), "\n";
}

print "2D Plot of #genes curves for each abundance vs ratio\n";

#plotting 3d surfaces
print "2D Plot data of #of events everywhere and within genes\n";
foreach my $rt (sort {$a <=> $b;} keys %ratios) {
    my %olapGenes = ();
#    print $rt, "-> ", scalar(@{$ratios{$rt}}), " ", scalar(@genes), "\n";
    my @rgenes = sort {&Gene::compareTo($a, $b);} (@{$ratios{$rt}});
    GenUtils::computeGeneOverlaps(\@genes, \@rgenes, \%olapGenes, 
				  0, 0, 0);
    
    print $rt*5, "\t", scalar(@{$ratios{$rt}}), " ", scalar(keys %olapGenes), "\n";
}

#plotting 3d surfaces
print "2D Plot data of #of events everywhere and within genes\n";
foreach my $rt (sort {$a <=> $b;} keys %ratios) {
    my %olapGenes = ();
#    print $rt, "-> ", scalar(@{$ratios{$rt}}), " ", scalar(@genes), "\n";
    my @rgenes = sort {&Gene::compareTo($a, $b);} (@{$ratios{$rt}});
    GenUtils::computeGeneOverlaps(\@genes, \@rgenes, \%olapGenes, 
				  0, 0, 0);
    
    print $rt*5, "\t", scalar(@{$ratios{$rt}}), " ", scalar(keys %olapGenes), "\n";
}


#plotting 3d surfaces
print "3D Plot data\n";
for(my $i = 0; $i < 4; $i++) {
    print "type ", $i, "\n";
    foreach my $abund (sort {$a <=> $b;} keys %{$abund[$i]}) {
	foreach my $rt (sort {$a <=> $b;} keys %{$ratios[$i]}) {
	    print $abund,"\t", $rt*5, "\t";
	    if(!defined($abund_ratios[$i]->{$abund}) || 
	       !defined($abund_ratios[$i]->{$abund}->{$rt})) {
		print "0\n";
	    }	
	    else {
		print $abund_ratios[$i]->{$abund}->{$rt}, "\n";
	    }
	}
    }
}

sub altSplScore {
    my($type, $abundance, $ratio) = @_;
    my $p_ab = 1;
    my $p_rt = 1;
    my $p_abrt = 1;
    
}
