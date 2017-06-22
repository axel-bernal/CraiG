#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $use_cid_offset = 0;
my $usage = "usage : $0 [--use-cid-offset] chunk_size < contigFile\n";

&GetOptions("-use-cid-offset" => \$use_cid_offset,
    );

my $chunkSize = shift || die "$usage\n";
my $c;
$/ = "\n>";

while(defined($c = Contig->new('File' => \*STDIN))) {
    my $base = 1;
    while($base < $c->{'Len'}) {
        my $id = $c->{'Id'};
        my $offset = 1;
        if($use_cid_offset && $id =~ /^(\S*[^\.]+)\.(\d+)$/) {
            $id = $1;
            $offset = $2;
        }

        my $substring = substr($c->{'Seq'}, $base - 1, $chunkSize);
        print ">$id.SPLIT_CONTIG.$c->{'Len'}.", $offset + $base - 1, "\n$substring\n";
        $base += $chunkSize;
    }
}

$/ = "\n";
