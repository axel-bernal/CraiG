#!/usr/bin/perl
#
# usage: genbankSplit.pl <data_set_name> <set #>
#
# this simple perl script splits up a GENBANK data set in $ARGV[0].xtrain.gbk/$ARGV[0].xtest.gbk sets to be used for x-evaluation procedures.
# $ARGV[0] contains data set WITHOUT extension, e.g. "combined_GB"
# the corresponding ".gbk" (or ".gbk.gz") and ".sets" files must exist.
#
# $ARGV[1] contains the number of the set to be used as the test set.
#
#

if ($#ARGV != 1) {
    die "Usage: $0 <data_set_name> <set #>";
}

open(LIST,"$ARGV[0].sets") || 
   die "Can't open list of data sets, $ARGV[0].sets";

$setnum = $ARGV[1];
$inset = 0;

while (<LIST>)
{
    chop;
    (/^$/) && ($inset=0);
    if ($inset)
    {
	s/ //g;
	$t{$_}=1;
    }
    if (/PART: $setnum/) 
    {
 	($inset=1);
    }
}

close(LIST);

open(TEST,">$ARGV[0].xtest.gbk");
open(TRAIN,">$ARGV[0].xtrain.gbk");
open(DATA,"< $ARGV[0].gbk") || open(DATA,"gunzip -c $ARGV[0].gbk.gz |") ||
    die "Can't open data file, $ARGV[0].gbk[.gz]";

$intest=0;
while (<DATA>)
{
    if (/^LOCUS/)
    {
	/^LOCUS\s+([^ ]+)/;
	if ($t{$1})
	{
	    $intest=1;
	}
	else {
	    $intest=0;
	}
    }
    if ($intest)
    {
	print TEST $_;
    }
    else
    {
	print TRAIN $_;
    }
}
	
