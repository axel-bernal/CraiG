#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Contig;

my $usage = "$0 maxLength < sequences\n";
my $maxLength = shift || die $usage;

my $c;
$/ = "\n>";
while(defined($c = Contig->new('File' => \*STDIN))) {
    my $base = 0;
    while($base < $c->{'Len'}) {
        print ">$c->{'Id'}.$base\n", substr($c->{'Seq'}, $base, $maxLength), "\n";
        $base += $maxLength;
    }
}
$/ = "\n";
