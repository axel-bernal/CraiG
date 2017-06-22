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

# Need to buffer the neg strand ouput as it's in reverse position
my @minus_buffer = ();
# keep track of last position in comp strand
my $last_comp_position = 0;

while(defined($line = <STDIN>)) {    
    if($line =~ />(\S+)/) { 
	if(scalar(@minus_buffer) > 0) {
	    print PREFIX_MINUS join("\n", @minus_buffer), "\n";
	    $last_comp_position = 0;
	}	    
	$contig = $1; 
	$curr_clen = $lengths{$contig};
	exit(1) if !defined($curr_clen);
	$factor = 1;
	@minus_buffer = ();
    } elsif($line =~ /\/\//) { 
	$factor = -1; 
    } elsif($line =~ /(\d+)\.\.(\d+)\s+(\S.+)$/) { 
	if($factor < 0) { 
	    if($1 > $last_comp_position) {
		unshift(@minus_buffer, "$contig\t".($curr_clen - $2)."\t".($curr_clen - $1 + 1)."\t".sprintf("%.3f", $factor*log($3)));
	    }
	    else {
		push(@minus_buffer, "$contig\t".($curr_clen - $2)."\t".($curr_clen - $1 + 1)."\t".sprintf("%.3f", $factor*log($3)));
	    }
	    $last_comp_position = $1;
	} else {
	    print PREFIX_PLUS "$contig\t", $1 - 1, "\t$2\t", sprintf("%.3f", $factor*log($3)), "\n";
	}
    } elsif($line =~ /(\d+)\s+(\S.+)$/) { 
	if($factor < 0) { 
	    if($1 > $last_comp_position) {
		unshift(@minus_buffer, "$contig\t".($curr_clen - $1)."\t".($curr_clen - $1 + 1)."\t".sprintf("%.3f", $factor*log($2)));
	    }
	    else {
		push(@minus_buffer, "$contig\t".($curr_clen - $1)."\t".($curr_clen - $1 + 1)."\t".sprintf("%.3f", $factor*log($2)));
	    }
	    $last_comp_position = $1;
	} else {
	    print PREFIX_PLUS "$contig\t", $1 - 1, "\t$1\t", sprintf("%.3f", $factor*log($2)), "\n";
	}
    }
}

if(scalar(@minus_buffer > 0)) {
    print PREFIX_MINUS join("\n", @minus_buffer), "\n";
}	    
									
close(PREFIX_PLUS);
close(PREFIX_MINUS);									     

