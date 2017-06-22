#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $usage = "usage : $0  < fastaFile";
my $c;

$/ = "\n>";
while(defined($c = Contig->new('File' => \*STDIN))) {
    if(!-d $c->{'Id'}) {
        mkdir($c->{'Id'});
    }
    else {
        next;
    }

    open(CONTIGS, ">$c->{'Id'}/contigs");
    open(NAME, ">$c->{'Id'}/name");
    print CONTIGS ">$c->{'Id'}\n$c->{'Seq'}\n";
    print NAME "$c->{'Legend'}\n";
    close(CONTIGS);
    close(NAME);
}
$/ = "\n";
