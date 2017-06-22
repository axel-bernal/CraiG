#!/usr/bin/perl -w
my $line;
while(defined($line = <STDIN>)) {
    my @line = split(/\s+/, $line);
#    print ">>>>",$line;
    if($line[2] eq "start_codon") {
	$line = <STDIN>;
	@line = split(/\s+/, $line);
	$line[2] = "initial-exon";
	$line[8] = "transgrp=".$line[11];
    }
    elsif($line[2] eq "stop_codon") {
	$forprint[-1]->[4] += 3;
	$forprint[-1]->[2] = ($#forprint == 0)? "single-exon":"final-exon";
#	print scalar(@forprint), " " , "a\n";
	for($i = 0; $i <= $#forprint; $i++) {
	    print join("\t", @{$forprint[$i]}), "\n";
	}
	@forprint = ();
	next;
    }
    else {
	$line[2] = "internal-exon";	
	$line[8] = "transgrp=".$line[11];
    }
#    print "***", join(" ", @line),"\n";
    push(@forprint, \@line);
}
