#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use GenUtils;
use Contig;

my $usage = "usage : $0 contigs < codr_xfasta > new_contigs\n";

my $contigs = shift ||  die $usage;
open(CONTIGS, "<$contigs");

my %codregs = %{&GenUtils::loadExtendedFasta(\*STDIN)};

$/ = "\n>";
while(defined($c = Contig->new('File' => \*CONTIGS))) {
    my $crArr = $codregs{$c->{'Id'}};
    my $len = 0;
    for(my $i = 0; $i < scalar(@$crArr); $i++) {
        my $cr = $crArr->[$i];
        if($cr->[0] == 0) {
            print ">$c->{'Id'}.$i\n",substr($c->{'Seq'}, $len, $cr->[1]), "\n";
        }
        $len += $cr->[1];
    }
}
$/ = "\n";

close(CONTIGS);
