#!/usr/local/bin/perl
#
# Author: David Kulp (dk/ucsc)
#
# usage: unpack.pl <data_set_name>
#
# this simple perl script creates a "unpacked" version of the data set
# in which a subdirectory is created, and each GenBank entry is a separate
# file.
#
# $ARGV[0] contains data set WITHOUT extension, e.g. "combined_GB"
#

if ($#ARGV != 0) {
    die "Usage: $0 <data_set_name>\n";
}

open(DATA,"< $ARGV[0].dat") || open(DATA,"gunzip -c $ARGV[0].dat.gz |") ||
    die "Can't open data file, $ARGV[0].dat[.gz]";

# why doesn't the permissions get set correctly???!!
mkdir("$ARGV[0]",0755);

while (<DATA>)
{
    if (/^LOCUS/)
    {
	/^LOCUS       ([^ ]+)/;
	close GENE;
	open(GENE,"> $ARGV[0]/$1");
    }
    print GENE $_;
}
	
