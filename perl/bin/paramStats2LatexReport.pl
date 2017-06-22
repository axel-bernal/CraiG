#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use strict;
use Getopt::Long;
my $usage = "$0 figure_name\n";
my $report = "both";
my $fig_name = shift || die "$usage";
my $line;
my $max_values = 0;
print "\\newpage\n\\begingroup\n\\scriptsize\n\\rowcolors{2}{gray!10}{}\\scriptsize\n\\begin{longtable}[!ht]{l|ll}\nFeature Name & Tag & Weight\\\\\n
\\hline\n";
$line = <STDIN>;
my @tbl_lines = ();
my $tbl_rows = 0;
while(defined($line = <STDIN>)) {
    my @line = split(/\s+/, $line);
    my @subscripts = ();
    my @superscripts = ();

    if($line =~ /\sValues$/) {
	$max_values = 1;
	print join("\n", @tbl_lines);
	print "\n\\hline\n\\caption{\\small Ordered Listing of Features and their Averaged Learnt Weights for \\".$fig_name."{} \\vspace*{.05in} \\newline\n\\scriptsize .}\n\\label{tbl:".$fig_name."_avgparamvals}\n\\end{longtable}\n\\endgroup\n";
	@tbl_lines = ();
	$tbl_rows = 0;
	print "\\newpage\n\\begingroup\n\\scriptsize\n\\rowcolors{2}{gray!10}{}\\scriptsize\n\\begin{longtable}[!ht]{l|ll}\nFeature Name & Tag & Weight\\\\\n";	
	next;
    }

    my $feat_name = " ";
    my $gram_len = -1;
    my $tag_start = 1;
    if($line[0] =~ /(\d)GRAM/) {
	$gram_len = $1;
	$feat_name = "count";
    }
    elsif($line[0] =~ /\-Signal/ || $line[0] =~ /\-PWM/ || $line =~ /QScore/) {
	if($line[0] !~ /QScore/) {
	    $feat_name = "PWM";
	}
	else {
	    my $sources = "";
	    $feat_name = "Lbins";
	    if($line[0] =~ /2Align/) {
		$sources = "{\\rm QScore}_{\\rm gap2}, {\\rm QScore}_{\\rm nraa}";
		$feat_name = "N\\text{-}".$feat_name;
	    }
	    elsif($line[0] =~ /^NR/ || $line[0] =~ /^St/) {
		$sources = "{\\rm QScore}_{\\rm nraa}";
	    }
	    else {
		$sources = "{\\rm QScore}_{\\rm gap2}";
	    }
	    $feat_name = "\\mathrm{".$feat_name."}_{25}({".$sources."})";		
	    if($line[0] =~ /\-Signal/ || $line[0] =~ /\-PWM/) {
		$feat_name = "{\\rm PWM}_{e} \\wedge ".$feat_name;
	    }
	}
    }
    else {
	print $line;
	next;
    }

    my $num_sources = -1;
    if($line[0] !~ /Union/ && $line[0] !~ /xPWM\-\dmer/ && $line[0] !~ /QScore/ && $line[0] =~ /(\d)mer/) {
	my $num_sources = $1;
	if($num_sources > 1) {
	    my $i = 0;
	    for(; $i < $num_sources; $i++) {
		my $source = $line[$i + 1];
		$source =~ s/\_/\\\_/g;
		push(@superscripts, "{\\rm ".$source."}");
	    }
	    $feat_name = $feat_name."cnj";
	    $tag_start = $i + 1;
	}
	else {
	    $tag_start = 2;
	    my $source = $line[1];
	    $source =~ s/\_/\\\_/g;
	    push(@subscripts, "{\\rm ".$source."}");
	}
    }

    if($gram_len >= 1) {
	push(@subscripts, $gram_len);
    }
    if($max_values) {
	push(@subscripts, "[".$line[-3]);
	push(@subscripts, $line[-2]."]");
    }

    if($line[0] =~ /\-Codonx/) {
	$feat_name = "d".$feat_name;
    }
    
    my $phase = -1;
    if($line[0] =~ /\-(\d)\S\S\-Frame$/) {
	$phase = $1;
	push(@subscripts, $phase);
    }

    if($line[0] =~ /PhaseHist/) {
	$phase = 0;
	push(@subscripts, $phase);
    }    

    if($line[0] !~ /QScore/) {
	$feat_name = "\\mathrm{".$feat_name."}";
    }
    for(my $i = 0; $i <= $#superscripts; $i++) {
	if($i == 0) {
	    $feat_name = $feat_name."^{";
	}
	$feat_name = $feat_name.$superscripts[$i];	
	if($i <= $#superscripts - 1) {
	    $feat_name = $feat_name.",";
	}
	elsif($i == $#superscripts) {
	    $feat_name = $feat_name."}";
	}
    }
    
    for(my $i = 0; $i <= $#subscripts; $i++) {
	if($i == 0) {
	    $feat_name = $feat_name."_{";
	}
	$feat_name = $feat_name.$subscripts[$i];	
	if($i <= $#subscripts - 1) {
	    $feat_name = $feat_name.",";
	}
	elsif($i == $#subscripts) {
	    $feat_name = $feat_name."}";
	}
    }
    
    if($line[0] =~ /Union/ || $line[0] =~ /xPWM\-\dmer/ && $line[0] !~ /QScore/) {
	$feat_name = $feat_name." \\ \\wedge \\ \\left(\\bigwedge_{e \\in \\mathcal E} {\\rm PWM}_{e}\\right)";
	$tag_start = $#line - ($max_values ? 4 : 2) + 1;
    }

    if($line[0] =~ /QScore/ && ($line[0] =~ /\-Signal/ || $line[0] =~ /\-PWM/)) {
	$feat_name = $feat_name.",\\ e \\in {\\mathcal E}";
    }
    if($line[0] =~ /Hist/) { # count duration feature
	$feat_name = "\\mathrm{bins}_w (\\log_2".$feat_name.")";
    }
    
    my $tbl_line = "\$".$feat_name."\$ & ";
    if($max_values) {
	$tag_start = $tag_start - $#line - 1;
	for(my $i = $tag_start; $i <= -4; $i++) { 
	    $tbl_line .= $line[$i]." ";
	}
    }
    else {
	$tag_start = $tag_start - $#line - 1;
	for(my $i = $tag_start; $i <= -2; $i++) { 
	    $tbl_line .= $line[$i]." ";
	}
    }	
    $tbl_line .= "& \$".(sprintf "%.2f", $line[-1])."\$ ";
    if($tbl_rows >= 40) {
	$tbl_lines[$tbl_rows - 40] .= " & ".$tbl_line."\\\\";
    }
    else {
	push(@tbl_lines, $tbl_line);
    }
    $tbl_rows++;	
}
print join("\n", @tbl_lines);
print "\n\\hline\n\\caption{\\small Ordered Listing of Individual Features and their Learnt Weights for \\".$fig_name."{} \\vspace*{.05in} \\newline\n\\scriptsize .}\n\\label{tbl:".$fig_name."_maxparamvals}\n\\end{longtable}\n\\endgroup\n";
