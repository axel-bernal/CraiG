#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use strict;

my $fwd = "";
my $bwd = "";

my $usage = "usage : $0 fwd_fasta bwd_fasta > all_fasta\n";

$fwd = shift || die $usage;
$bwd = shift || die $usage;
open(FWD, $fwd) || die $usage;
open(BWD, $bwd) || die $usage;

my($line1, $line2);
$/ = "\n>";
while(defined($line1 = <FWD>)) {
#    print "\"", $line1,"\"";
    $line1 =~ /^>?(\S+)\s+([\n\s\S]+)\n>?$/;
    my $id1 = $1;
    my $seq1 = $2;

    $line2 = <BWD>;
    defined($line2) || die $id1;
    $line2 =~ /^>?(\S+)\s+\d+\n([\n\s\S]+)\n>?$/;
    ($id1 eq $1) || die $id1;
    print ">$id1\t$seq1\n";
    print "//\n$2\n";
}
$/ = "\n";

close(FWD); close(BWD);
