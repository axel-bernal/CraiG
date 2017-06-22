#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use strict;

my $i;
my $id;
my $usage = "$0 prefix comma-sep-indexes < signalP\noutputs SPgff format for indexes specificied in the argument list\n";

my $sigprefix = shift || die $usage;
my $indexes = shift || die $usage;
my %idsf = ();
my %idsb = ();
my @line = ();

my @indexes = split(/\,/, $indexes);

while(<STDIN>) {
    @line = split(/\s+/);
    my @details = split(/\,/, $line[8]);
    foreach my $ind (@indexes) {
#	print $details[$ind], "\n";
	if($details[$ind] =~ /([^\=]+)\=([0-9,\.,\-]+)$/) {
	    open(SIGPOUT, ">>$sigprefix.$1.SPgff") || die "Can't open $sigprefix.$1 for output";
	    print SIGPOUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$2\t",
	    "$line[6]\t$line[7]\t$details[0]\n";
	    close(SIGPOUT);
	}
    }
}
