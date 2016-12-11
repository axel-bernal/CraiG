#!/usr/bin/perl -w                                                                          
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-                                                                              

use lib "$ENV{CRAIG_HOME}/perl/lib";
use enum qw(STRAND_FWD STRAND_COMP);
use GenUtils;
use strict;

my $usage = "$0 < list_junction_files > table\n";
my $line;
while(defined($line = <STDIN>)) {
    next if $line =~ /^\#/;
    chomp $line;
    $line =~ /(\S+)\s+(\S+)\s+(\S+)/;
    my $file = $1;
    my $dataset = $2;
    my $sample = $3;
    open(JUNCT, $file) || die "Can't open $file\n";
    my $line = <JUNCT>;
    while(1) {
	my $junction = GenUtils::readLocsFromString($line, undef);
	if(!defined($junction)) {
	    last;
	}

	my $exons = $junction->{'Exons'};
	next if $#$exons <= 0;
	print "$dataset\t$sample\t$exons->[0]->[0]\t";
#	print join("\t", @{$exons->[0]}),"\t\t", join("\t", @{$exons->[1]}),"\n";
	if($junction->{'Strand'} == STRAND_FWD) {
	    print $exons->[0]->[2] + 1, "\t", $exons->[1]->[1] - 1;
	}
	else {
	    print $exons->[1]->[1] + 1, "\t", $exons->[0]->[2] - 1;
	}
	print "\t$exons->[0]->[3]\t$exons->[0]->[3]\t0\t";
	print $junction->{'Strand'} == STRAND_FWD ? "+":"-", "\n";
        $line = <JUNCT>;
    }
    close(JUNCT);
}
