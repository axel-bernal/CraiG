#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use GenUtils;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < fasta > xfasta\n";
my $c;

$/ = "\n>";
while(defined($c = Contig->new('File' => \*STDIN))) {
    my $seq = &GenUtils::fasta2xfasta($c->{'Id'}, $c->{'Seq'});
    print ">$c->{'Id'}\n$seq";
}
$/ = "\n";
