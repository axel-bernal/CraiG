#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Getopt::Long;
use strict;

my $usage = "usage : $0 < sigPosFile > outSigPosFile\n";

my $offset = 50;
my $line=  <STDIN>;
while(defined($line) && $line =~ /^>(\S+)\s+(\d+)/) {
    my $id = $1;
    my $length = $2;
    while(defined($line = <STDIN>) && $line !~ /^>/) {
        if($line =~ /^2\s+(\d+)/) {
            print "D ", $1 + $offset,"\n";
        }
        elsif($line =~ /^3\s+(\d+)/) {
            print "A ", $1 + $offset,"\n";
        }
    }
    $offset += 50 + $length;
}
