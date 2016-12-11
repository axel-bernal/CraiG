#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use strict;

$| = 1;
my $usage = "usage: $0 [-partial|-all] < wormMart > gtf_file\n";
my $partial = 0;
my $all = 0;
&GetOptions("-partial" => \$partial,
            "-all" => \$all
            );

my $line;
my $transcript_id;
my %parents = ();
my $offset;
my %confirmed = ();
my $transcript;
my %genes = ();
my @lines = ();

while(defined($line = <STDIN>)) {
    if($line =~ /^\#.*FASTA/) {
	last;
    }

    if($line =~ /^\#/) {
        next;
    }
    my @line = split(/\t/, $line);
    
    if($line[2] =~ /gene/i && $line[8] =~ /ID=(Gene\:)?([^\;,\s]+)/) {
        $genes{$2} = 1;
    }
    if($line[2] =~ /^mRNA$/ || $line[2] =~ /^transcript$/) {
        $line[8] =~ /ID=(Transcript\:)?([^\;,\s]+)/;        
        $transcript = $2;
        $line[8] =~ /Parent=(Gene\:)?([^\;,\s]+)/;
        my $gene = $2;

        if(defined($genes{$gene})) {
            $parents{$transcript} = $gene;
        }
        else { warn "Gene $gene was not defined in $line[8]\n"; }

        if($partial && $line[8] =~ /status=Partially_confirmed/) {
            $confirmed{$transcript} = 1;
        }
        elsif($line[8] =~ /status=Confirmed/) {
            $confirmed{$transcript} = 1;
        }

    }
    
    if($line[2] !~ /CDS/i && $line[2] !~ /five\_prime\_UTR/i && $line[2] !~ /three\_prime\_UTR/i && $line[2] !~ /exon/ && $line[2] !~ /start_codon/ && $line[2] !~ /stop_codon/) {
        next;
    }

    if($line[2] !~ /UTR/ && $line[7] eq ".") {
	$line[7] = 0;
    }
    push(@lines, \@line);
}

foreach $line (@lines) {
    my @transcripts = ();
    my $detail = $line->[8];
    my $exon_number = -1;
    if($detail =~ /Parent=(Transcript)?/) {
        while($detail =~ /Parent=(Transcript\:)?([^\;,\s]+)\,?(.*)$/) {
            push(@transcripts, $2);
            $detail = "Parent=$3";
        }
    }
    elsif($detail =~ /ID=(CDS\:)?([^\;,\s]+).*$/) {
        push(@transcripts, $2);
    }    

    if($line->[2] =~ /exon/) {
	$line->[8] =~ /ID=[^\-]+\-([0-9]+)\;/;
	$exon_number = $1;
    }
    foreach $transcript (@transcripts) {
        if(!defined($confirmed{$transcript}) && !$all) {
            next;
        }

#        $offset = 0;
#    if($line[3] > 1000000) {
#        $offset = int($line[3]/1000000)*1000000;
#        $line[3] += -$offset;
#        $line[4] += -$offset;
#    }
        
        if(!defined($parents{$transcript})) {
            $parents{$transcript} = $transcript;
        }
#    print $line[0], ":", $offset + 1, "..", $offset + 1000000, "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5], "\t", $line[6], "\t", $line[7], "\tgene_id \"", $parents{$transcript}, "\"; transcript_id \"", $transcript, "\";\n";
        print $line->[0], "\t", $line->[1], "\t", $line->[2], "\t", $line->[3], "\t", $line->[4], "\t", $line->[5], "\t", $line->[6], "\t", $line->[7], "\tgene_id \"", $parents{$transcript}, "\"; transcript_id \"", $transcript;
	if($exon_number >= 0) {
	    print "\"; exon_number \"", $exon_number;
	}
	print "\";\n";
    }
}
