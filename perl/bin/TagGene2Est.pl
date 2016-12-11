#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
$file = shift ||  die "$0 in_file";
`sed s/SINGLE_EXON/EST_EXON/g $file.ge > $file.1`;
`sed s/LAST_EXON/EST_EXON/g $file.1 > $file.2`;
`rm $file.1`;
`sed s/INIT_EXON/EST_EXON/g $file.2 > $file.3`;
`rm $file.2`;
`sed s/LONG_INTRON/INTRON/g $file.3 > $file.4`;
`rm $file.3`;
`sed s/INTRON/EST_INTRON/g $file.4 > $file.est.ge`;
`rm $file.4`;
