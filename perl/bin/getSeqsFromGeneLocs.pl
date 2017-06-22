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

#open(AA, "a");
#my @list = ();
#while(<AA>) {
#    my $line = $_;
#    chomp($line);
#    push(@list, $line);
##

#foreach my $id (keys %contigs) {
#    $contigs{$id}->{'Seq'} =~ s/\s+//g;
#    print ">$id\n", $contigs{$id}->{'Seq'}, "\n";
#}
#__END__;
my %genes = % { $org->loadGenes(\*STDIN)};
foreach my $id (keys %genes) {
    $genes{$id}->dump(\*STDOUT);
}
