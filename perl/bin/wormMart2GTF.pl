#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Getopt::Long;

$| = 1;
my $usage = "usage: $0 [-CDS|-PARTIAL|-COMPLETE] < wormMart > gtf_file\n";

my $line;
my $len = 0;
my $transcript_id;
my $offset;
my $coord = 10;
my $cds = 0;
my $partial = 0;
my $confirmed = 0;

&GetOptions("-CONFIRMED" => \$confirmed,
            "-PARTIAL" => \$partial,
            "-CDS" => \$cds,
	    );

if($cds) {
    $coord = 12;
}

while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    chomp($line);
    
    my @line = split(/\t/, $line);

    if($partial && $line[5] !~ /Partially_confirmed/ && $line[5] !~ /Confirmed/) {
        next;
    }
    elsif($confirmed && $line[5] !~ /Confirmed/) {
        next;
    }

    $line[8] =~ /(\S+)\.exon(\d+)$/;
    $transcript_id = $1;

    if($2 == 1) {
        $len = 0;
    }

    if($line[$coord] eq "" || $line[$coord + 1] eq "") {
        next;
    }

    $offset = 0;
    if($line[$coord] > 1000000) {
#        $offset = int($line[7]/1000000)*1000000;
#        $line[7] += -$offset + 1;
#        $line[8] += -$offset + 1;
    }
    
    print $line[6];
#    print ":", $offset + 1, "..",$offset + 1000000;
    print "\tWormMart\tCDS\t", $line[$coord],"\t", $line[$coord + 1], "\t.\t", ($line[7] == 1 ? "+" : "-"), "\t", $len % 3, "\t", "gene_id \"", $line[1], "\"; transcript_id \"", $transcript_id, "\";\n";
    $len += abs($line[$coord] - $line[$coord]) + 1;
}
