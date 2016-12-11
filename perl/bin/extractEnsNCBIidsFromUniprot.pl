#!/usr/bin/perl
use strict;
my $line;
while(defined($line = <STDIN>)) {
    my $id;
    $line =~ /[^\:]+\:(\S+)\s+(\S+)\s+/;
    my $swiss_id = $1;
    my $gene_id = $2;
    
    if($line =~ /RefSeq\s+\S+\s+(NM\_\S+)\s+/) {
	$id = $1;
    }
    elsif($line =~ /Ensembl\s+(ENST\S+)\;\s+/) {
	$id = $1;
    }
    else {
	warn "$line does not have either Ensemble nor RefSeq ids\n";
	last; 
    }
    print "$id\t$swiss_id\t$gene_id\n";
}

