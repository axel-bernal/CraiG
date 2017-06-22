#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use strict;

my $threshold = 0.99;
my $usage = "usage : $0 -t threshold[0 .. 0.99] < blastFile > new_genes";

&GetOptions("-t=f" => \$threshold,
            );

my %matches = ();
my $line;
select STDOUT; $| = 1;

$/ = "\nBLASTN";
while(defined($line = <STDIN>)) {
    chomp $line;
    $line =~ /\nQuery=\s*(\S+)\s*\n\s*\((\d+)\s+letters[.,\d,\S,\s]+\nSequences[^\n]+\n\n([^>]+)(>[.,\d,\S,\s,\n]+)$/;

    my $id = $1;
    my $length = $2;
    print "$id $length\n";
    defined($id && $length && $3 && $4) || next;
    my $matches = &extractBlastPMatches($3, $4, $threshold);
    
    for(my $i = 0; $i < scalar(@{$matches}); $i++) {
        my $match = $matches->[$i];
        $match->[3] || next;
        print "debug \"$match->[0]\" \n\"$match->[1]\"", $match->[2]*1.0/$match->[3], "\n";
        
        if($match->[2]*1.0/$length < $threshold || 
           defined($matches{$id})) {
            next;
        }
        $matches{$id} = $match
    }
}

$/ = "\n";

foreach my $gene (keys %matches) {
    if(defined($matches{$gene})) { 
        print ">$gene ", join(" ", @{$matches{$gene}}),"\n";
    }
}

sub extractBlastPMatches {
    my($blastPString, $blastIdent, $score) = @_;
    my $remnant = $blastPString;
    my @results = ();
    my @remnant = split(/\n/, $remnant);
    my $i = 0;

    while($i < scalar(@remnant)) {
        if($remnant[$i] !~ /\S/) {
            $i++;
            next;
        }
        if($remnant[$i] =~ /(\S+)\s+\d+\s+([^\s]+)\s+\d+$/) {
            print "debug match -> $1 $2\n";
            push(@results, [$1, $2, 0.0, 0.0]);
            $i++;
        }
        else {
            last;
        }
    }

    $remnant = $blastIdent;
    @remnant = split(/>/, $remnant);
    shift @remnant;
    $i = 0;
    while($i < scalar(@remnant)) {
        print "remnant \"", $remnant[$i],"\"\n";
        
        if($remnant[$i] =~ /^(\S+)[.,\s,\S,\d]+\n\sIdentities\s+\=\s+(\d+)\/(\d+)\s*/) {
            my $id = $1;
            $results[$i]->[2] = $2;
            $results[$i]->[3] = $3;
            print "remn $id $2 $3 $results[$i]->[0]\n";
            ($id eq $results[$i]->[0]) || die "$id $results[$i]->[0]\n";
            $i++;
        }
        else {
            last;
        }
    }
      
    return \@results;
}

