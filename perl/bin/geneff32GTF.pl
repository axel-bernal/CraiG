#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
$| = 1;
my $usage = "usage: $0 < wormMart > gtf_file\n";

my $line;
my $transcript_id;
my %parents = ();
my $offset;

while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    my @line = split(/\t/, $line);
#    print "****", $line[8], "\n";
    if($line[8] =~ /ID=Transcript\:(\S+);Parent=Gene\:(\S+)\;/) {
        $parents{$1} = $2;
    }

    if($line[2] !~ /CDS/i) {
        next;
    }
    $line[8] =~ /Parent=Transcript\:(\S+)\;/;
    $offset = 0;
    if($line[3] > 1000000) {
        $offset = int($line[3]/1000000)*1000000;
        $line[3] += -$offset;
        $line[4] += -$offset;
    }

    print $line[0], ":", $offset + 1, "..", $offset + 1000000, "\t", $line[1], "\t", $line[2], "\t", $line[3], "\t", $line[4], "\t", $line[5], "\t", $line[6], "\t", $line[7], "\tgene_id \"", $parents{$1}, "\"; transcript_id \"", $1, "\";\n";
}
