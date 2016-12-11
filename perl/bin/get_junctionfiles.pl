#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use strict;
my $top_level="/eupath/data/EuPathDB/workflows/ToxoDB/19/data/";
my $usage = "$0 list_of_RNA_Seq_libraries genome_assembly_vers\n";
my $list = shift || die $usage;
my $version = shift || die $usage;
open(LIST, "$list") || die "Can't open $list\n";
my $line;

while(defined($line = <LIST>)) {
    $line =~ /(\S+)\s+(\S+)/;
    my $dbname = $1;
    my $sample = $2;
    $dbname =~ /([^\_]+)\_(\S+)/;
    my $org = $1;
    my $dataset = $2;
    my $junct_file = "/gpfs/fs121/h/abernal/Data/toxodb-19/AltSplice/me49-$version/$dataset"."_"."$sample.preproc/tgondii-rna.$dataset.$sample.chr.junction.locs\t$dataset\t$sample";
    print $junct_file,"\n";
}
