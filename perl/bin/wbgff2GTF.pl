#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use strict;

$| = 1;
my $usage = "usage: $0 [-include-partial|-include-predicted] < wormbasegff > gtf_file\n";
my $incl_partial = 0;
my $incl_pred = 0;
&GetOptions("-include-partial" => \$incl_partial,
            "-include-predicted" => \$incl_pred,
            );

my $line;
my $transcript_id;
my %parents = ();
my $offset;
my %included = ();
my $transcript;
my %genes = ();

while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    my @line = split(/\t/, $line);
    
    if($line[2] =~ /Transcript/i) {
        $line[8] =~ /Transcript\s\"([^\"]+)\"/;
        my $transcript = $1;
        $line[8] =~ /Gene\s\"([^\"]+)\"/;
        my $gene = $1;
        $genes{$gene} = 1;
        $parents{$transcript} = $gene;
        if($incl_pred || 
           ($incl_partial && $line[8] =~ /confirmed\"/i) ||
           ($line[8] =~ /\"Confirmed\"/)) {
            $included{$transcript} = 1;
        }
    }
    
    if($line[2] !~ /coding_exon/i && $line[2] !~ /five\_prime\_UTR/i && $line[2] !~ /three\_prime\_UTR/i) {
        next;
    }

#    print "$line[2] $line[8]";

    my @transcripts = ();
    $line[8] =~ /Transcript\s\"([^\"]+)\"/;
    push(@transcripts, $1);
    
    foreach $transcript (@transcripts) {

        if(!defined($included{$transcript})) {
            next;
        }
        
#        $offset = 0;
#    if($line[3] > 1000000) {
#        $offset = int($line[3]/1000000)*1000000;
#        $line[3] += -$offset;
#        $line[4] += -$offset;
#    }
        
        if(!defined($parents{$transcript})) {
            warn "$transcript has no parents\n";
            $parents{$transcript} = $transcript;
        }
#    print $line[0], ":", $offset + 1, "..", $offset + 1000000, "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5], "\t", $line[6], "\t", $line[7], "\tgene_id \"", $parents{$transcript}, "\"; transcript_id \"", $transcript, "\";\n";
        print $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5], "\t", $line[6], "\t", $line[7], "\tgene_id \"", $parents{$transcript}, "\"; transcript_id \"", $transcript, "\";\n";
    }
}
