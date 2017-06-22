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

my $org = Organism->new();
my %marks = %{&GenUtils::loadExtendedFasta(\*STDIN)};
my %done = ();

foreach my $id (keys %marks) {
    my $id2=  $id;
    my $strand = undef;
    if($id =~ /([^\+,\-,\s]+)(\+|\-)/) {
        $id2 = $1;
        $strand = $2;
    }
    if($done{$id2})  {
        next;
    }
    if(defined($strand)) {
        print ">$id2\n";
        foreach my $t (@{$marks{$id2."+"}}) {
            ($t->[1] >= 0) or die $id;
            print "$t->[0] x $t->[1]\n";
        }
        print "\/\/\n";
        foreach my $t (@{$marks{$id2."-"}}) {
            ($t->[1] >= 0) or die $id;
            print "$t->[0] x $t->[1]\n";
        }
    }
    else {       
        print join("\n", @{$marks{$id2}}), "\n";
    }
    $done{$id2} = 1;
}


