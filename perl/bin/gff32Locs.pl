#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use strict;
use File::Temp;

my $contigFile;
my $exon_tag = "";
my $cds_tag = "";
my $sort = 0;
my $preserve_phases = 0;
my $usage = "usage : $0 [--preserve-phases] -cds-tag CDS_TAG -exon-tag EXON_TAG -ctg contigFile < gff3_file\nThis script only works if you redirect the gff3_file directly to STDIN\n";

&GetOptions(
    "--preserve-phases" => \$preserve_phases,
    "-sort" => \$sort,
    "-cds-tag=s" => \$cds_tag,
    "-exon-tag=s" => \$exon_tag,
    "-ctg=s" => \$contigFile,
    );

my (undef, $filename1) = File::Temp::tempfile();
my (undef, $filename2) = File::Temp::tempfile();

die "$usage\nTag for cds or exon must be given" if $cds_tag eq "" && $exon_tag eq "";
die "$usage\nContigs must be provided\n" if $contigFile eq "";

if($cds_tag ne "") {
#    open local(*STDIN), '<', "whatever" or die $!;
    system("cat - | grep $cds_tag | gff32GTF.pl -all | gtf2Locs.pl | fixTranslatedProduct.pl ".($preserve_phases ? "--preserve-phases " : "")."-ctg $contigFile | grep \">\" | sort > $filename1");
#    close(STDIN) or die $!;
}

if($exon_tag ne "") {
    seek(STDIN, 0, 0);
    system("cat - | grep $exon_tag | gff32GTF.pl -all | gtf2Locs.pl | sort > $filename2");
}

if($cds_tag ne "" && $exon_tag ne "") {
    my $sort_opt = ($sort ? "-sort" : "");
    print `produceUTRannot.pl $sort_opt -utr $filename2 < $filename1`;
}
elsif($cds_tag ne "") {
    print `cat $filename1`;
}
elsif($exon_tag ne "") {
    print `cat $filename2`;
}
