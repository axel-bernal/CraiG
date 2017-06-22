#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-my %lines = ();

my $line;
while(defined($line = <>)) {
    if(defined($lines{$line})) {
    }
    else {
        $lines{$line} = 1;
        print $line;
    }
}
