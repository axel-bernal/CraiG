#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use strict;
$| = 1;
my $usage = "usage: $0 < wormMart > gtf_file\n";
my $line;
my $transcript_id;
my %genes = ();
my %transcripts = ();
my %children = ();
my $offset;
my $g;
my $t;
my $id = 1;
my $i;

print "##gff-version 3\n";
while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    my @line = split(/\s+/, $line);
    if($line[2] !~ /CDS/ && $line[2] !~ /UTR/ && 
       $line[2] !~ /start_codon/ && $line[2] !~ /stop_codon/) {
        next;
    }
    $line[9] =~ /\"(\S+)\"/;
    $g = $1;
    $line[11] =~ /\"(\S+)\"/;
    $t = $1;
    
    my @new_line = ();
    for($i = 0; $i < 9; $i++) {
        push(@new_line, $line[$i]);
    }       
    
    if($new_line[0] =~ /([^\:]+)\:(\d+)\.\.\d+/) {
	$new_line[0] = $1;
	$new_line[3] += ($2 - 1);
	$new_line[4] += ($2 - 1);
    }

    $new_line[8] = "ID=$id;Parent=Transcript:$t;Name=$t";

    if(!defined($genes{$g})) {
        my @gene_line = ();
        push(@gene_line, @new_line);
        $gene_line[2] = "gene";
        $gene_line[5] = $new_line[5];
        $gene_line[7] = '.';
        $gene_line[8] = "ID=Gene:$g;Name=$g";
        $genes{$g} = \@gene_line;
    }
    else {
        if($genes{$g}->[3] > $new_line[3]) {
            $genes{$g}->[3] = $new_line[3];
        }
        if($genes{$g}->[4] < $new_line[4]) {
            $genes{$g}->[4] = $new_line[4];
        }
    }
    if(!defined($transcripts{$t})) {
        @{$transcripts{$t}} = ();
        if(!defined($children{$g})) {
            @{$children{$g}} = ();
        }
        push(@{$children{$g}}, $t);
        my @trans_line = ();
        push(@trans_line, @new_line);
        $trans_line[2] = "mRNA";
        $trans_line[5] = $new_line[5];
        $trans_line[7] = ".";
        $trans_line[8] = "ID=Transcript:$t;Parent=Gene:$g;Name=$t";
	if(scalar(@line) == 13) {
	    $trans_line[8] .= ";".$line[12];
	}
        push(@{$transcripts{$t}}, \@trans_line);
    }
    else {
        if($transcripts{$t}->[0]->[3] > $new_line[3]) {
            $transcripts{$t}->[0]->[3] = $new_line[3];
        }
        if($transcripts{$t}->[0]->[4] < $new_line[4]) {
            $transcripts{$t}->[0]->[4] = $new_line[4];
        }
	if($transcripts{$t}->[0]->[5] ne ".") {
	    $transcripts{$t}->[0]->[5] += $new_line[5];
	}
    }

    push(@{$transcripts{$t}}, \@new_line);
    $id++;
}

foreach my $g (keys %children) {
    if($genes{$g}->[5] ne ".")  {
	$genes{$g}->[5] = 0;
	foreach $t (@{$children{$g}}) {
	    for($i = 1; $i < scalar( @{$transcripts{$t}}); $i++) {
		$genes{$g}->[5] += $transcripts{$t}->[$i]->[5];
	    }
	}
    }

    print join("\t", @{$genes{$g}}), "\n";
    
    foreach $t (@{$children{$g}}) {
        print join("\t", @{$transcripts{$t}->[0]}), "\n";
        for($i = 1; $i < scalar( @{$transcripts{$t}}); $i++) {
            print join("\t", @{$transcripts{$t}->[$i]}), "\n";
        }
    }
}
