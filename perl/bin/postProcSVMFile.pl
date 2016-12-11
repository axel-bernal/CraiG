#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Genefinder;
use Contig;
use strict;

my $usage = "usage :$0 fastaFile scoreFile <  svmFile\n";
my $fastaFile = shift || die "$usage";
my $scoreFile = shift || die "$usage";
my $line = "";
open(SCORE, "<$scoreFile");
open(FASTA, "<$fastaFile");

my @linesSVM = ();
while($line = <STDIN>) {
    chomp($line);
    my @line = split(/\s+/, $line);
    push(@linesSVM, \@line);
}

my %linesSC = ();
while($line = <SCORE>) {
    chomp($line);
    my @line = split(/\s+/, $line);
    $linesSC{$line[0]} = \@line;
}

my @linesFAS = ();
while($line = <FASTA>) {
    if($line !~ /^>/) {
	next;
    }
    chomp($line);
    my @line = split(/\s+/, $line);
    push(@linesFAS, \@line);
}

for (my $i = 0; $i < scalar(@linesSVM); $i++) {
    print join(" ", @{$linesSVM[$i]}), " ";
    my $num = scalar(@{$linesSVM[$i]});
    if(!defined($linesSC{$linesFAS[$i]->[0]})) {
	print STDERR $linesFAS[$i]->[0],"\n";
	next;
    }
    print $num, ":", $linesSC{$linesFAS[$i]->[0]}->[1], " ",
    $num + 1, ":", $linesSC{$linesFAS[$i]->[0]}->[2], " ", 
    $num + 2, ":", $linesSC{$linesFAS[$i]->[0]}->[3], " ",
    $num + 3, ":", $linesSC{$linesFAS[$i]->[0]}->[4], "\n";
}

close(FASTA);
close(SCORE);
