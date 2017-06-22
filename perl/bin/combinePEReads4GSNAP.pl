#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

my $usage = "$0 file1 file2 > gsnap_output\n";

my $file1 = shift || die $usage;
my $file2 = shift || die $usage;
open(FILE1, "<$file1") || die $usage;
open(FILE2, "<$file2") || die $usage;
my $line1;
my $line2;
while(defined($line1 = <FILE1>)) {
    $line2 = <FILE2>;
    if($line1 =~ /^\@(.+)\/\d/) {
	$line1 = <FILE1>;
	$line2 = <FILE2>;
	print ">$1\n$line1$line2";
    }
}
