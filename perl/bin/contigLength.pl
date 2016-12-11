#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile = "";
my $id = "";
my $subContigs = 0;
my $usage = "usage : $0 -l -ctg contigFile -id contigId\n";

&GetOptions("-id=s" => \$id,
	    "-ctg=s" => \$contigFile,
            "-l" => \$subContigs,
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

open(CONTIGS, "<$contigFile") || die "file $contigFile not found\n";
my $c;
my $length;
if($subContigs) {
    my %contigs = %{&GenUtils::readLocsFormat(\*CONTIGS, 0, 'CONTIGS', undef)};
    if($id eq "") {
        foreach $c (keys %contigs) {
            $length = $contigs{$c}->{'Exons'}->[0]->[2] -
                $contigs{$c}->{'Exons'}->[0]->[1] + 1;
            print $length, "\t$c\n";
        }
    }
    else {
        $length = $contigs{$id}->{'Exons'}->[0]->[2] - 
            $contigs{$id}->{'Exons'}->[0]->[1] + 1;
        print $length, "\n";
    }        
}

if($id eq "") {
    $/ = "\n>";
    while(defined($c = Contig->new('File' => \*CONTIGS))) {
        $length = length($c->{'Seq'});
        print $length, "\t$c->{'Id'}\n";
    }
    $/ = "\n";
}
else {
    my $org = Organism->new();
    my %contigs = % {$org->loadContigs(\*CONTIGS)};
    $length = length($contigs{$id}->{'Seq'});
    print $length, "\n";
}

close(CONTIGS);
