#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < contigFile";
my $c;

$/ = "\n>";
while(defined($c = Contig->new('File' => \*STDIN))) {
    my $offset = 1;
    my $new_cid = $c->{'Id'};
    if($c->{'Id'} =~ /^([^\.]+)\.(\d+)$/) {
        $new_cid = $1;
        $offset = $2;
    }
    
    print ">$c->{'Id'}\t$new_cid", "_", "$offset", "_", length($c->{'Seq'}) + $offset - 1, "\n";
}
$/ = "\n";
