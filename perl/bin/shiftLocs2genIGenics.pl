#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use GenUtils;
use Contig;

my $usage = "usage : $0 codr_xfasta < orig_locs\n";
my $codr = shift || die $usage;
open(CODR, "$codr") || die "Can't open $codr\n";
my %codregs = %{&GenUtils::loadExtendedFasta(\*CODR)};
close(CODR);

my %genes = % {&GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};

foreach my $key (keys %genes) {
    my $gene = $genes{$key};
    my $crArr = $codregs{$gene->{'Exons'}->[0]->[0]};;
    my $len = 0;
    my $offset = 0;
    for(my $i = 0; $i < scalar(@$crArr); $i++) {
        my $cr = $crArr->[$i];
        if($cr->[0] == 1) {
            $offset = $len;
        }
        $len += $cr->[1];
    }
    ($len != 0) || die $gene->{'Id'};
#    print STDERR "$gene->{'Id'} $offset $len\n";
    $gene->shiftCoords($offset, $len);
    $gene->dump(\*STDOUT, undef, 1);
}
