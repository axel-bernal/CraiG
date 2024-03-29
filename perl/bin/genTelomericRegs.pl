#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $contigFile;
my $tel_len = 0;
my $usage = "usage : $0 --tel-len TEL_LENGTH(Kb) -ctg contigFile > out_contig_file\n";

&GetOptions("-tel-len=i" => \$tel_len,
            "-ctg=s" => \$contigFile,
    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my $tlen = $tel_len*1000;
my $tel_seq = "TTAGGG" x ($tlen / 6);
foreach my $id (keys %contigs) {
    my $clen = length($contigs{$id}->{'Seq'});
    my $cseq = $tel_seq.substr($contigs{$id}->{'Seq'}, $tlen, $clen - $tlen).$tel_seq;
    print ">$id\n$cseq\n";
}

