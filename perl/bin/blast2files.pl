#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use Getopt::Long;
use strict;

my $contigFile;
my $ds = 0;
my $ups = 0;
my $dir = "";
my $usage = "usage : $0 -dir directory < fastaFile";

&GetOptions("-dir=s" => \$dir,
	    );

(length($dir) != 0) || die $usage;
my $line;

if(! -d "$dir") {
    mkdir("$dir");
}

$/ = "\nBLASTN";
while(defined($line = <STDIN>)) {
    chomp $line;
    if($line !~ /^BLASTN/) {
        $line = "BLASTN".$line;
    }
    $line =~ /\nQuery=\s*(\S+)\s*\n/;
    my $id = $1;
    if(! -e "$dir/$id.blast") {
        `touch "$dir/$id.blast"`;
    }
    open(CONTIGS, ">>$dir/$id.blast");
    print CONTIGS "$line\n";
    close(CONTIGS);
}
$/ = "\n";
