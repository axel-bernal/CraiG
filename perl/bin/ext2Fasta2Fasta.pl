#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);
my $ual = 0;
my $usage = "usage : $0 [-u] < ext2fasta > fasta\nOptions\n\t-u\tfo including unaligned regions\n";

&GetOptions(
            "-u" => \$ual,
            );

my $line =  <STDIN>;
my $fillChar;

while(defined($line) && $line  =~ />(\S+)\s+(\d+)\s+(\S)/) {
    my $cid = $1;
    my $clen = $2;
    $fillChar = $3;
    my $seq = $fillChar x $clen;
    print ">$cid\n";

    while(defined($line = <STDIN>) && $line  !~ />(\S+)\s+(\d+)\s+(\S)/) {
        if($line =~ /^(\d+)\s+(\S+)$/) {
            my $pos = $1;
            my $this_seq = $2;
            my $len = length($this_seq);

	    for(my $i = 0; $i < $len; $i++) {
		substr($seq, $pos + $i, 1) = substr($this_seq, $i, 1);
	    }
        }
        elsif($line =~ /\/\//) {
            die "Only forward strand filters allowed\n";
        }
    }

    print $seq,"\n";
}
