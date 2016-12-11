#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Contig;
use Getopt::Long;
use strict;

my $dir = "";
my $prefix = "";
my $append = 0;
my $usage = "usage : $0 [-append] [-dir directory] [-prefix prefix] < fastaFile";

&GetOptions("-dir=s" => \$dir,
	    "-append" => \$append,
	    "-prefix=s" =>\$prefix,
	    );

die $usage if $dir eq "" && $prefix eq "";
my $fn = ($dir ne "" ? "$dir/$0" :"$prefix.$0");

if($dir ne "") {
    unlink("$dir") unless $append;
    mkdir("$dir") unless (-e $dir);
}

open(OC, ">$fn");
while (<STDIN>)
{
    my $i = 0; 
    if(/^>(\S+)\s*/) {
	close(OC);
	$fn = ($dir ne "" ? "$dir/$1" : "$prefix.$1");
	$fn = ">$fn" unless !$append;
	open(OC, ">$fn") || die "Can't write in $fn\n";
    }
    print OC $_;
}
close(OC);
unlink($dir ne "" ? "$dir/$0" : "$prefix.$0");
