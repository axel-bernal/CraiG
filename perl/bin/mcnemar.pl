#!/usr/bin/perl -w
my $usage = "$0 ann_file pred_file other_pred_file [-exon\-transcript]\n";
my $annot = shift || die "$usage\n";
my $pred = shift || die "$usage\n";
my $opred = shift || die "$usage\n";
my $type = shift || die "$usage\n";

if($type eq "exon") {
    if(! -f "$annot.exons") {
	system("$ENV{CRAIG_HOME}/perl/bin/listExons.pl < $annot.locs | sort > $annot.exons");    
    }
    if(! -f "$pred.exons") {
	system("$ENV{CRAIG_HOME}/perl/bin/listExons.pl < $pred.locs | sort > $pred.exons");
	system("join $pred.exons $annot.exons > $pred.correx");
	system("join -v 1 $pred.exons $annot.exons > $pred.incorrex");
    }
    if(! -f "$opred.exons") {
	system("$ENV{CRAIG_HOME}/perl/bin/listExons.pl < $opred.locs | sort > $opred.exons");
	system("join $opred.exons $annot.exons > $opred.correx");
	system("join -v 1 $opred.exons $annot.exons > $opred.incorrex");
    }
    
    print `join $pred.correx $opred.correx | wc`;
    print `join -v 1 $pred.correx $opred.correx | wc`;
    print `join -v 2 $pred.correx $opred.correx | wc`;
    print `join $pred.incorrex $opred.incorrex | wc`;
}
else {
    system("biointers.py $pred.locs $annot.locs -by ll > $pred.corrt");
    system("biodiff.py $pred.locs $annot.locs -by ll > $pred.incorrt");
    system("biointers.py $opred.locs $annot.locs -by ll > $opred.corrt");
    system("biodiff.py $opred.locs $annot.locs -by ll > $opred.incorrt");
    
    print "c = ", `biointers.py $pred.corrt $opred.corrt -by ll | wc`;
    print "b = ", `biodiff.py $pred.corrt $opred.corrt -by ll | wc`;
    print "d = ", `biodiff.py $opred.corrt $pred.corrt -by ll | wc`;
    print "a = ", `cat $pred.incorrt | filterOlapGenes.pl -o $opred.locs | wc`;
    print "e = ", `cat $opred.incorrt | filterOlapGenes.pl -o $pred.locs | wc`;
#    print "a = ", `filterOlapGenes.pl -o $annot.locs < $opred.locs | filterOlapGenes.pl -o $pred.locs | wc`;
#    print "e = ", `filterOlapGenes.pl -o $annot.locs < $pred.locs | filterOlapGenes.pl -o $opred.locs | wc`;
    print `biointers.py $pred.incorrt $opred.incorrt -by ll | wc`;
    print `biodiff.py $pred.incorrt $opred.locs -by ll| wc`;
    print `biodiff.py $opred.incorrt $pred.locs -by ll| wc`;      
    print "add b+e for 01 and d+a for 10\n";
}
