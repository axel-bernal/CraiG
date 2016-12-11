#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use strict;

my $threshold = 1e-10;
my $usage = "usage : $0 -t threshold < blastFile > new_genes";

&GetOptions("-t=f" => \$threshold,
            );

my %discarded = ();
my $line;
select STDOUT; $| = 1;

$/ = "\nBLASTP";
while(defined($line = <STDIN>)) {
    chomp $line;
    $line =~ /\nQuery=\s*(\S+)\s*[.,\d,\S,\s]+\nSequences[^\n]+\n\n([^>]+)/;
#    print "debug \"$1\" \n\"$2\" \n";
    my $id = $1;
    defined($id && $2) || die $id;
    my $matches = &extractBlastPMatches($2, $threshold);
    my $selfMatch = 0;
    $discarded{$id} = "no";
    
    foreach my $match (@{$matches}) {
        if($match->[0] eq $id) {
            $selfMatch = 1;
            if($discarded{$id} eq "yes") {
                last;
            }
            next;
        }

        if($match->[1] > 1.0e-20 || 
           (defined($discarded{$match->[0]}) && $discarded{$match->[0]} eq "yes")) {
            next;
        }
        $discarded{$id} = "yes";
        if($selfMatch) {
            last;
        }
    }
    
    if(!$selfMatch) {
        warn "Can't find selfMatch for $id, seg has masked it and it will be discarded\n";
        $discarded{$id} = "yes";
    }
}

$/ = "\n";

foreach my $gene (keys %discarded) {
    if($discarded{$gene} eq "no") {
        print ">$gene\n";
    }
}

sub extractBlastPMatches {
    my($blastPString, $score) = @_;
    my $remnant = $blastPString;
    my @results = ();
    my $match = "";

    while(1) {
        if($remnant =~ /^([^\n]+)\n([.,\d,\S,\s]*)$/) {
            $match = $1;
            $remnant = $2;
        }
        else {
            last;
        }
        
        if($match =~ /(\S+)\s+\S+\s+\d+\s+([^\s]+)\s+\d+$/) {
#            print "debug $match -> $1 $2\n";
            push(@results, [$1, $2]);
        }
        else {
            die $match;
        }
    }
    
    return \@results;
}

