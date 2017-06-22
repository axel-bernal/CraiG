#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use GenUtils;
use Utils;
use strict;
use enum qw(STRAND_FWD STRAND_COMP);

my $clens_file = "";
my $filter_is_score = 0;
my $cols_are_pos  = "";
my $contig_file = ""
my $usage = "usage : $0 [--is-score-filter --scores-are-positions=col1[,col2..coln] --ctg subcontiglocs --clen contig_lens < subcontig_sigfilter > contig_sigfilter. Input should be ordered byt contig id\n";

&GetOptions("-ctg=s" => \$contig_file,
            "-scores-are-positions=s" => \$cols_are_pos,
            "-is-score-filter" => \$filter_is_score,
	    "-contig-len=s" => \$clen_file;
    );

die $usage unless $clens_file ne "" > 0 && $contig_file ne "";

open(SUBCONTIGS, "<$contig_file") || die "file $contig_file not found\n";
my %subcontigs = % { &GenUtils::readLocsFormat(\*SUBCONTIGS, 0, 'TRAIN', undef)};
close(SUBCONTIGS);

my %contig_lens = ();
open(LENGTHS, "$clensFile") ||  die $usage;
my $line;
while(defined($line = <LENGTHS>)) {
    $line =~ /(\d+)\s+(\S+)/;
    $contig_lens{$2} = $1;
}
close(LENGTHS);

my %subfilterVals = ();
my %filterVals = ();

my $def_val = ($filter_is_score ? 0 : '-');

$line = <STDIN>;
while(1) {
    my $subcid = GenUtils::loadFilterVals4Id(\*STDIN, \$line, \%subfilterVals, $def_val);
    if($subcid eq "") {
	last;
    }

    (defined($subcontigs{$subcid})) || die "$subcid does not exist\n";
    my $contig = $contigs->{$subcid}->{'Exons'}->[0];
    my $cid = $contig->[0];
    my $clen = $contig_lens{$cid};
    GenUtils::createFilterValArrays(\%filterVals, $cid, $clen, $def_val);

    my $subclen = $contig->[2] - $contig->[1] + 1;
    for(my $strand = 0; $strand < 2; $strand++) {
	my $fvals = $subfilterVals{$subcid}->[$strand];
	for(my $i = 1; $i <= $subclen; $i++) {
	    my @vals = split(/\s+/, $fvals->[$i]);

	    if($filter_is_score) {
		
	    }
	    else {
		my $pos = $contig->[1] + $i - 1;
		$filterVals{$cid}->[$strand]->[$pos] = 
		    for(my $k = 1; $k <= $#fvals; $k += 3) {
			$max = $fvals[$k+2] unless $max >= $fvals[$k+2];
		}
	    }
	}
    }
}

GenUtils::printScoreFilterVals(\%filterVals, \*STDOUT, STRAND_COMP);
print STDERR "Done\n";
