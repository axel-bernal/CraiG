#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
my $usage = "$0 lengths prefix\n";

my $lengths = shift || die "No lengths file given\n".$usage;
my $prefix = shift || die "No prefix given\n".$usage;
my $line;

open(LENGTHS, $lengths) || die "Can't open file $lengths";

while(defined($line = <LENGTHS>)) {
    $line =~ /^(\S+)\s+(\d+)$/;
    $lengths{$1} = $2;
}
close(LENGTHS);

my $contig;
my $factor;
my $curr_clen;
open(PREFIX_PLUS, ">$prefix"."_plus.bed") || die "Can't open file $prefix"."_plus.bed\n";
open(PREFIX_MINUS, ">$prefix"."_minus.bed") || die "Can't open file $prefix"."_minus.bed\n";

while(defined($line = <STDIN>)) {    
    if($line =~ />(\S+)/) { 
	$contig = $1; 
	$curr_clen = $lengths{$contig};
	exit(1) if !defined($curr_clen);
	$factor = 1;
    } elsif($line =~ /\/\//) { 
	$factor = -1; 
    } elsif($line =~ /(\d+)\.\.(\d+)\s+(\S.+)$/) { 
	if($factor < 0) { 
	    print PREFIX_MINUS "$contig\t", $curr_clen - $2, "\t", $curr_clen - $1 + 1, "\t", $factor*$3, "\n";
	} else {
	    print PREFIX_PLUS "$contig\t", $1 - 1, "\t$2\t", $factor*$3, "\n";
	}
    } elsif($line =~ /(\d+)\s+(\S.+)$/) { 
	if($factor < 0) { 
	    print PREFIX_MINUS "$contig\t", $curr_clen - $1, "\t", $curr_clen - $1 + 1, "\t", $factor*$2, "\n";
	} else { 
	    print PREFIX_PLUS "$contig\t", $1 - 1, "\t$1\t", $factor*$2, "\n";
	} 
    }
}
									
close(PREFIX_PLUS);
close(PREFIX_MINUS);									     

