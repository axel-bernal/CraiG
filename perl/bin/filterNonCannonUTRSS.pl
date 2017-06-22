#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

my $usage = "$0 locs contigs utrs\n";
my $locs = shift || die $usage;
my $contigs = shift || die $usage;
my $utrs = shift || die $usage;
`sort $locs > $locs.$$`;
`filterNonCannonSS.pl -ctg $contigs -filter gene < $utrs > $utrs.$$ 2> $utrs.ssreport`;
`grep ">" $utrs.$$ | grep "5UTR" | sed 's/.5UTR//'  | cut -f1 | sort > a.$$`;
`grep ">" $utrs.$$ | grep "3UTR" | sed 's/.3UTR//'  | cut -f1 | sort > b.$$`;
`join a.$$ b.$$ | sort > c.$$`;
`join -v 1 a.$$ c.$$ > aa.$$`;
`mv aa.$$ a.$$`;
`join -v 1 b.$$ c.$$ > aa.$$`;
`mv aa.$$ b.$$`;
`join $locs.$$ a.$$ | strip3UTR.pl | sort > a.locs.$$`;
`join $locs.$$ b.$$ | strip5UTR.pl | sort > b.locs.$$`;
`join $locs.$$ c.$$ | sort > c.locs.$$`;
`cat a.locs.$$ b.locs.$$ c.locs.$$ | sort > d.locs.$$`;
`grep -v "\\[" $locs | grep -v "\\]" > e.locs.$$`;
print `cat d.locs.$$ e.locs.$$ | filterNonCannonSS.pl -ctg $contigs -filter gene 2> $locs.ssreport | grep ">" |  filterRepeatedGenes.pl`;
`rm $locs.$$ $utrs.$$ a.$$ b.$$ c.$$ a.locs.$$ b.locs.$$ c.locs.$$ d.locs.$$ e.locs.$$`;

