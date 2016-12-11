#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Getopt::Long;
use Contig;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

$| = 1;
my $usage = "usage: $0 -ctg contigs [-est] < alignment_file > locs_file\n";

my $contigFile = "";
my $est = 0;
my $checkStart = 0;
my $canonical = 0;
my $filterStops = 0;

&GetOptions("-est" => \$est,
	    "-check-start" => \$checkStart,
	    "-canonical" => \$canonical,
	    "-filter-inframe-stops" => \$filterStops,
            "-ctg=s" => \$contigFile);

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
$/ = "\n>";
my $c = Contig->new('File' => \*CONTIG);
close(CONTIG);
$/ = "\n";    
my $line;
my $counter = 0;

while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    chomp($line);
    
    my @line = split(/\s/, $line);
    my $id = $c->{'Id'}.".".$line[4];
    my $hasStart = 0;
    my $hasStop = 0;
    my $hasAcceptor = 0;
    my $hasDonor = 0;
    my $strand = STRAND_FWD;

    if($line[1] < $line[0]) {
        $strand = STRAND_COMP;
    }

    my $acceptor = $c->getSequence(2, $line[0], $line[0], 0, 0, $strand);
    if(length($acceptor) == 3 && ($acceptor =~ /^ag/i)) {# || $acceptor =~ /^ac/i)) {
        $hasAcceptor = 1;
    }

    my $donor = $c->getSequence(0, $line[1], $line[1], 2, 0, $strand);
    if(length($donor) == 3 && ($donor =~ /gt$/i)) {# || $donor =~ /gc$/i || $donor =~ /at$/i))  {
        $hasDonor = 1;
    }

    my $acceptStart = 0;
    my $acceptStop = 0;
    my $phase = 0;
    my $alseq = $c->getSequence(0, $line[0], $line[1], 3, 0, $strand);
    my $numStops = length($alseq);

    if(!$est) {
	my $starts = GenUtils::phaseOccurrences($alseq, $phase, "atg");
	my $stopPhase = -1;
	my @stops = ();
#	print STDERR $id,".",$counter+1," ", $#{$starts}, " ", length($alseq), "\n";
	for($phase = 0; $phase < 3; $phase++) {
	    push(@stops, GenUtils::phaseOccurrences($alseq, $phase, "taa", "tga", "tag"));

	    if(scalar(@{$stops[$phase]}) > 0 && ${$stops[$phase]}[-1] == length($alseq) - 3) {
		$hasStop = 1;
		$stopPhase = $phase;
	    }
	}

	if($hasStop && scalar(@{$stops[$stopPhase]}) == 1) {
	    $acceptStop = 1;
	}

#	print STDERR "stop $acceptStop $stopPhase\n";
	if($#$starts >= 0 && $starts->[0] == 0) {
	    $hasStart = 1;
	    if(scalar(@{$stops[0]}) == 0 || ($hasStop == 1 && $stopPhase == 0 && scalar(@{$stops[0]}) == 1)) {
		$acceptStart = 1;
		if($checkStart) {
		    if(($line[2] == 1 && $line[3] > $line[2]) || 
		       ($line[3] == 1 && $line[2] > $line[3])) {
		    }
		    else {
			print STDERR "Wrong coordinates for $id match.. perhaps not a start?\n";
		    }
		}

		if($acceptStop && $stopPhase != 0) {
		    if($hasAcceptor) {
			$acceptStop = 0;
		    }
		    else {
			print STDERR "start and stop have different phase\n";
		    }
		}
	    }
	}

	if($acceptStop && ${$stops[$stopPhase]}[-1] != length($alseq) - 3) { 
	    if($strand == STRAND_FWD) {
		$line[1] -= 3;
	    }
	    else {
		$line[1] += 3;
	    }
	}

	if(!$acceptStart && !$acceptStop) { #no translation signals, find feasible phase
	    $numStops = length($alseq);
	    my $minPhase = 0;
	    for($phase = 0; $phase < 3; $phase++) {
		if(scalar(@{$stops[$phase]}) < $numStops) {
		    $minPhase = $phase;
		    $numStops = scalar(@{$stops[$phase]});
		}
	    }
	    $phase = $minPhase;

	    if($numStops) {
		print STDERR "\tphase $phase gives $numStops stops for $id.",$counter+1,"\n";
		if($filterStops) {
		    print STDERR "removing\n";
		    next;
		}
	    }	    
	}
    }
    else {
    }

    my $loc = $c->{'Id'}."_".$line[0]."_".$line[1];

    if(!$acceptStart) {
	if($hasAcceptor || !$canonical) {
	    $loc = "<".$loc;
	    if($phase) {
		$loc = $phase.$loc;
	    }	    
	}
	else {
	    print STDERR "no acceptor for $id.",$counter+1,"\n";
	    next;
	}
    }

    if(!$acceptStop) {
	if($hasDonor || !$canonical) {
	    $loc = $loc.">";
	}
	else {
	    print STDERR "no donor for $id.",$counter+1,"\n";
	    next;
	}
    }

    $counter++;
    print ">$id.$counter\t$loc\t$line[5]\t$line[6]\n"; 
}
