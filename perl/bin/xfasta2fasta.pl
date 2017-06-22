#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  

my $usage = "usage : $0 < xfasta > fasta\n";

$_ = <>;
while(defined($_) && $_ =~ /^>/) {
    print $_;
    my $seq = "";
    while(defined($_ = <>) && $_ !~ /^>/) {
        $_ =~ /^([^\s,>])\s+x\s+(\d+)/;
        if($2) {
            $seq .= $1 x $2;
        }
    }
    if($seq ne "" ) {
        print $seq,"\n";
    }
}
