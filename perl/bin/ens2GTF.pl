#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Getopt::Long;

$| = 1;
my $usage = "usage: $0 < ens_file > gtf_file\n";

my $line;
my $cdslen = 0;
my $transcript_id = "";
my $cds = 0;
my $utr = 0;
my @last_line = ();

while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    chomp($line);
    
    my @line = split(/\t/, $line);
    
    $cds = 1;
    $utr = 0;
    if($line[3] eq "" || $last_line[9] ne $line[9]) {
        if($last_line[9] ne $line[9]) {
            ($line[-2] == $line[1]) || warn "different transcript init for $line[9] -> $line[-2] $line[1]\n";
            ($last_line[-1] == $last_line[2]) || warn "different transcript end for $last_line[9] -> $last_line[-1] $last_line[2]\n";
        }
        $cdslen = 0;
    }

    @last_line = @line;

    if($line[3] eq "") {
        $cds = 0;
        $utr = 1;
    }
    else {
        if($line[3] > $line[1]) {
            $line[2] = $line[3] - 1;
            $utr = 1;
        }

        if($line[6] < $line[2]) {
            $line[1] = $line[6] + 1;
            $utr = 1;
        }
    }

#    $offset = 0;
    if($line[1] > 1000000) {
#        $offset = int($line[7]/1000000)*1000000;
#        $line[7] += -$offset + 1;
#        $line[8] += -$offset + 1;
    }
    
#    print ":", $offset + 1, "..",$offset + 1000000;

    if($utr) {
        print "chr".$line[10]."\tVegaKnown\t";
        print "UTR\t", $line[1],"\t", $line[2], "\t.\t", ($line[11] == 1 ? "+" : "-"), "\t0\t", "gene_id \"", $line[8], "\"; transcript_id \"", $line[9], "\";\n";
    }
    if($cds) {
        print "chr".$line[10]."\t";
        print "\tVegaKnown\tCDS\t", $line[3],"\t", $line[6], "\t.\t", ($line[11] == 1 ? "+" : "-"), "\t", $cdslen % 3, "\t", "gene_id \"", $line[8], "\"; transcript_id \"", $line[9], "\";\n";
        $cdslen += abs($line[6] - $line[3]) + 1;
    }
}
