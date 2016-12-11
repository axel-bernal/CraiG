#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

my $maxNum = 200000;
$maxNum = shift || die;
my $sum = 0;
for(my $i = 1; $i < $maxNum; $i++) {
    $sum = $sum + 1/$i;
    if(!($i % 100)) {
        print $i," ", $sum, " ", log($i), "\n";
    }
}
