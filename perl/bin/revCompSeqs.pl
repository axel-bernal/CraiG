#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";  
use GenUtils;

# -*- perl -*-

# usage:  revCompSeqs.pl  < dna_seq  > reversed_dna_seq
#

my $header;
my $id;
my $line = <>;
while(defined($line) && $line =~ /^>(\S+)(\s?.+)$/) {
    $id = $1;
    $header = $2;
    my @lines = ();
    while(defined($line = <>) && $line !~ /^>.+/) {
        push(@lines, $line);
    }
    my $seq = join("",@lines);
    $seq =~ s/\s+//g;
    print ">comp",$id, $header, "\n", &GenUtils::revComp($seq),"\n";
}
