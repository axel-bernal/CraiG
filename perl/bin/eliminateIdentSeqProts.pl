#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < proteinSeqs_w_Header\n";

my $org = Organism->new();

my %contigs = % {$org->loadContigs(\*STDIN)};
my %prots = ();

foreach my $id (keys %contigs) {
    $contigs{$id}->{'Seq'} =~ s/\s+//g;
#    if(!defined($prots{$contigs{$id}->{'Seq'}}) || $id =~ /annotation/) {
        $prots{$contigs{$id}->{'Seq'}} = $id;
#    }
}

foreach my $prot (keys %prots) {
    print ">$prots{$prot}\n";
}
