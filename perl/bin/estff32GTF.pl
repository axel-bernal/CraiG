#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Utils;
use strict;
$| = 1;
my $usage = "usage: $0 < ests_in_gff3_fmt > gtf_file\n";

my $line;
my $tid;


while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    my @line = split(/\t/, $line);

    if($line[5] < 90) {
        next;
    }
    $tid = undef;
    if($line[8] =~ /Target=Sequence\:([^\;,\s]+)/) {
        $tid = $1;
    }
    if(!defined($tid)) {
        next;
    }
    print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t\"gene_id \"", $tid, "\"; transcript_id \"", $tid, "\";\n";

}
