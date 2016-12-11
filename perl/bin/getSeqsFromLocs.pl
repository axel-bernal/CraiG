#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $ds = 0;
my $ups = 0;
my $usage = "usage : $0 -ups upstream -ds downstream -ctg contigFile < locations\n";

&GetOptions("-ups=i" => \$ups,
            "-ds=i" => \$ds,
            "-ctg=s" => \$contigFile,
            );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my $line;
while(defined($line = <STDIN>) && ($line =~ /^>\S+\s+(.+)\s+?(\d+)?b?p?$/)) {
    my @locs = @{&GenUtils::readLocExons($1)};
    print $line;
    foreach my $loc (@locs) {
        print substr($contigs{$loc->[0]}->{'Seq'}, $loc->[1] - 1, $loc->[2] - $loc->[1] + 1);
    }
    print "\n";
}
