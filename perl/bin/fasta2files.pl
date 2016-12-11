#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $dir = "";
my $n = 1;
my $usage = "usage : $0 [-n num_seqs[1]] -dir directory < fastaFile";

&GetOptions("-dir=s" => \$dir,
            "-n=i" => \$n,
	    );

(length($dir) != 0) || die $usage;
my $c;
mkdir("$dir");
$/ = "\n>";
$c = Contig->new('File' => \*STDIN);

while(defined($c)) {
    my $i = 0; 
    open(CONTIGS, ">$dir/$c->{'Id'}");
    while($i < $n) {
        print CONTIGS ">$c->{'Id'}\t$c->{'Legend'}\n$c->{'Seq'}\n";
        $c = Contig->new('File' => \*STDIN);
        if(!defined($c))  {
            last;
        }
        $i++;
    }
    close(CONTIGS);    
}
$/ = "\n";
