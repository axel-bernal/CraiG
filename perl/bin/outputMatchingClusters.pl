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
my $use_scores = 0;
my $report_gtf = 0;
my $preserve_ids = 0;
my $usage = "usage : $0 [--report-gtf] [--use-scores] [--match-cds annotation] < clusters\n";

&GetOptions("--match-cds=s" => \$annotation,
	    "--use-scores" => \$use_scores,
	    "--report-gtf" => \$report_gtf,
	    "--preserve-ids" => \$preserve_ids,
    );

my $base_annot = ($annotation eq "" ? "" : basename($annotation));
my $line = <STDIN>;
my %counters = ();

while(defined($line) && $line =~ /^\/\/\s+(\S+)\s+(\S.+)$/) {
    my @annotations = split(/\s+/, $2);
    my $header_line = $line;
    my $genes = GenUtils::readLocsFormat2Array(\*STDIN, 0, 'TRAIN', \$line);
#    print "-> \"", $base_annot, "\" ", join(" ", @annotations), " ", scalar(@$genes), "\n";

    if(($base_annot eq "" && scalar(@annotations)) || 
       ($base_annot ne "" && scalar(@annotations) > 1 && $annotations[0] eq $base_annot)) {	
	my $header_detail = "";
	my $max5 = -1;
	my $max3 = -1;
	my $max = -1;
	my $max5_score = -10000000;
	my $max3_score = -10000000;
	my $max_score = -10000000;
	my $max5_length = -1;
	my $max3_length = -1;
	my $max_length = -1;
	my $offset = ($base_annot ne "" ? 1 : 0);
	
	for(my $i = $offset; $i < scalar(@$genes); $i++) {
	    my $gene = $genes->[$i]->[1];
	    my $m5utrs =  $gene->{'5UTRExons'};
	    my $m3utrs =  $gene->{'3UTRExons'};
	    my $cdsExons = $gene->{'Exons'};
	    my $length5 = 0;
	    my $length3 = 0;
	    my $lengthCDS = 0;
	    my $tPscore = scalar(@$m3utrs);
	    my $fPscore = scalar(@$m5utrs);
	    my $CDSscore = scalar(@$cdsExons);

	    if($use_scores) {
		$tPscore = $tPscore ? $gene->getThreePWeight() : -10000000;
		$fPscore = $fPscore ? $gene->getFivePWeight() : -10000000;
		$CDSscore = $CDSscore ? $gene->getCDSWeight() : -10000000;
	    }

	    if($base_annot eq "") { 
		if($i > 0) {
		    $header_detail .= ";";
		}

		$header_detail .= "Sample=".$annotations[$i];		
		$header_detail .= "\;Score=".($fPscore + $CDSscore + $tPscore);
		if(scalar(@$m5utrs)) {
		    $header_detail .= "\;FiveUTR_Score=".$fPscore;
		    $header_detail .= "\;FiveUTR=".scalar(@$m5utrs);
		}
	    }
	    
	    for(my $j = 0; $j < scalar(@$m5utrs); $j++) {
		my $utr = $m5utrs->[$j];
		$length5 += abs($utr->[2] - $utr->[1]) + 1;
		if($base_annot eq "" && $j == 0) {
		    if($gene->{'Strand'} == STRAND_FWD) { 
			$header_detail .= $utr->[1]."..".$utr->[2];
		    }
		    else {
			$header_detail .= $utr->[2]."..".$utr->[1];
		    }
		}
	    }

	    if(updateMax(\$max5_score, \$max5_length, 
			 $fPscore, $length5) || $i == $offset) {
		$max5 = $i;
	    }

	    for(my $j = 0; $j < scalar(@$cdsExons); $j++) {
		my $exon = $cdsExons->[$j];
		$lengthCDS += abs($exon->[2] - $exon->[1]) + 1;
	    }
		
	    if($base_annot eq "" && scalar(@$m3utrs)) {
		$header_detail .= "\;ThreeUTR_Score=".$tPscore;
		$header_detail .= "\;ThreeUTR=".scalar(@$m3utrs);
	    }		    

	    for(my $j = 0; $j < scalar(@$m3utrs); $j++) {
		my $utr = $m3utrs->[$j];
		$length3 += abs($utr->[2] - $utr->[1]) + 1;
		if($base_annot eq "" && $j == 0) {
		    if($gene->{'Strand'} == STRAND_FWD) { 
			$header_detail .= $utr->[1]."..".$utr->[2];
		    }
		    else {
			$header_detail .= $utr->[2]."..".$utr->[1];
		    }
		}
	    }
	    
	    if(updateMax(\$max3_score, \$max3_length, 
			 $tPscore, $length3) || $i == $offset) {
		$max3 = $i;
	    }
	    my $this_length = $length5 + $lengthCDS + $length3;
	    my $this_score = $tPscore + $CDSscore + $fPscore;

	    if(updateMax(\$max_score, \$max_length,
			 $this_score, $this_length) || $i == $offset) {
		$max = $i;
	    }	    

#	    print "DEBUG $i $gene->{'Id'} $this_score $max $max3_score $max5_score $this_length\n";
	}

	if($max >= 0) {
	    if($base_annot ne "") {
		if(scalar(@{$genes->[$max5]->[1]->{'5UTRExons'}})) {
		    $header_detail = "FiveUTR_Sample=$annotations[$max5];FiveUTR_Score=$max5_score";
		}
		if(scalar(@{$genes->[$max3]->[1]->{'3UTRExons'}})) {
		    if(length($header_detail)) {
			$header_detail .= ";";
		    }
		    $header_detail .= "ThreeUTR_Sample=$annotations[$max3];ThreeUTR_Score=$max3_score";
		}
		print "// ", $header_detail,"\n";
		$genes->[0]->[1]->dump(\*STDOUT, undef, 1);
		$genes->[$max]->[1]->{'Id'} = $genes->[0]->[1]->{'Id'}.".CR".$max;
		$genes->[$max]->[1]->{'3UTRExons'} = $genes->[$max3]->[1]->{'3UTRExons'};
		$genes->[$max]->[1]->{'5UTRExons'} = $genes->[$max5]->[1]->{'5UTRExons'};
	    }
	    else {
		print "// ", $header_detail,"\n";
		$genes->[$max]->[1]->{'Id'} = $annotations[$max].".gene".$counters{$annotations[$max]}++ if !$preserve_ids;
	    }

	    if($report_gtf) {
		&GenUtils::writeGeneGTFFormat($genes->[$max]->[1], \*STDOUT, $header_detail);
	    }
	    $genes->[$max]->[1]->dump(\*STDOUT, undef, 1);
	}
    }
}

sub updateMax {
    my ($m_score, $m_length, $score, $length) = @_;
	
    if($score  > $$m_score) {
	$$m_score = $score;
	$$m_length = $length;
	return 1;
    }
    elsif($score == $$m_score) {
#	print "DEBUG2 $score $$m_score $length $$m_length\n";
	if($length > $$m_length) {
	    $$m_length = $length;
	    return 1;
	}
    }
    return 0;
}
