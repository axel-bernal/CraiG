#!/usr/bin/perl
#
# usage: split.pl <type> <data_set_name> <set #>
#
# this simple perl script splits up a type data set in train/test sets.
# $ARGV[0] contains data set 
# $ARGV[1] contains the number of the set to be used as the test set.
#
# writes data out to "data_set_name.test" and "data_set_name.train"
#

($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";

if ($#ARGV != 3) {
    die "Usage: $0 <fasta|gtf|genbank> file sets <set #>.\n";
}

my %types = ('fasta' => '^>(\S+)\s?',
             'gtf' => '^(\S+)\.gb\.fa',
             'genbank' => '^LOCUS\s+([^ ]+)',
             );

open(LIST,"$ARGV[2]") || 
    die "Can't open list of data sets, $ARGV[2]";

$setnum = $ARGV[3];
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
    if (/PART: $setnum$/) 
    {
	($inset=1);
    }
}

close(LIST);

open(DATA,"<$ARGV[1]") || die "Can't open data file, $ARGV[0]";

$intest=0;
while (<DATA>)
{
    if(/$types{$ARGV[0]}/) {
	if ($t{$1}) {
	    $intest=1;
	}
	else {
	    $intest=0;
	}
    }
    if ($intest)
    {
	print $_;
    }
}
	
