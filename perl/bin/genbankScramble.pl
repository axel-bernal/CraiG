#!/usr/bin/perl
#
# this programs takes a GSDB-type GenBank flatfile and outputs a set of
# test/train sets.
#
# usage: genbankScramble.pl <GenBankfile>

if ($#ARGV!=0) { die "usage: scramble.pl <GenBankfile>";}

$NUMPARTS = 8;

@IDS = `grep "LOCUS" $ARGV[0]`;

$size = $#IDS;

srand;

for ($i=0; $i<=$size; $i++)
{
    $swap = int(rand($size+1));
    $tmp=$IDS[$i];
    $IDS[$i] = $IDS[$swap];
    $IDS[$swap] = $tmp;
}

$i=0;
for ($part=0; $part<$NUMPARTS; $part++)
{
    print "\nPART: $part\n";
    for ($j=0; $j<($size/$NUMPARTS); $j++,$i++)
    {
	$IDS[$i] =~ /^LOCUS\s+([^ ]+)/;
	print "$1\n";
    }
}

#put possible remaining few in final bin
for (;$i<=$size;$i++)
{
    $IDS[$i] =~ /^LOCUS\s+([^ ]+)/;
    print "$1";
}

    
