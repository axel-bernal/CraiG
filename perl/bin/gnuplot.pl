#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
                     
$x = "[]";
$y = "[-2:4]";                                                 
$usage = "usage: $0 id xrange yrange\n";
$id = shift || die $usage;
if(defined($ARGV[0])) {
    $x = shift;
}
if(defined($ARGV[0])) {
    $y = shift;
}
open(TMP,">$id.plot") || die "Cannot open $id.plot\n";
#print TMP "set size 2, 1\n";
print TMP "set ylabel \'Feature Scores\'\n";
print TMP "set xlabel \'bps\'\n";
#print TMP "set output \'$id.png\'\n";
print TMP "plot $x $y '$id' u 1:2 t \"start\" w impulses, '$id' u 1:3 t \"stop\" w impulses, '$id' u 1:4 t \"donor\" w impulses, '$id' u 1:5 t \"acceptor\" w impulses, '$id' u 1:6 t \"fr1\" w l, '$id' u 1:7 t \"fr2\" w l, '$id' u 1:8 t \"fr3\" w l, '$id' u 1:9 t \"real\" w l 4\n";
close TMP;
`gnuplot -persist $id.plot`;
