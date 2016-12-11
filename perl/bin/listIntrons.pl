#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $usage = "usage : $0 --include-ids < locations\n";
my $include_ids = 0;
&GetOptions("-include-ids" => \$include_ids);
    
my %genes = % {&GenUtils::readLocsFormat(\*STDIN)};
foreach my $id (keys %genes) {
    listIntrons($include_ids, $genes{$id}->{'3UTRExons'}, 
		$genes{$id}->{'Strand'}, \*STDOUT);
    listIntrons($include_ids, $genes{$id}->{'Exons'}, 
		$genes{$id}->{'Strand'}, \*STDOUT);
    listIntrons($include_ids, $genes{$id}->{'5UTRExons'},
		$genes{$id}->{'Strand'}, \*STDOUT);
}
 
sub listIntrons {
    my($include_ids, $exons, $strand, $file) = @_;
    for(my $i = 0; $i < scalar(@{$exons}); $i++) {
        if($i == 0) {
            next;
        }

	my $ctg = $exons->[$i]->[0];
	my $ibeg = $exons->[$i - 1]->[2] + 1;
	my $iend = $exons->[$i]->[1] - 1;
        if($strand == STRAND_COMP) {
	    $ibeg = $exons->[$i - 1]->[2] - 1;
	    $iend = $exons->[$i]->[1] + 1;
	}

	if($include_ids) {
	    print ">", $ctg, ":", $ibeg, "..", $iend, "\t";
	}
	
	print $ctg, "_", $ibeg, "_", $iend;

	if(defined($exons->[$i]->[3])) {
	    print "\t", $exons->[$i]->[3];
	}
	print "\n";
    }
}
