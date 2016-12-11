#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  

my $usage = "usage :$0  < listIds\n";
my $line = "";
my @gc = (43, 51, 57, 100);
my @ids = ();
my @types = ("coding", "intron");

while(defined($line = <STDIN>)) {
    chomp($line);
    push(@ids, $line);
}

#delete
`cp genes.fsa a`;
foreach my $id (@ids) {
    `grep -w -v $id a > a0`;
    `cp a0 a`;
}
`mv a0 genes.fsa`;
`cp eval_genes.fsa a`;
foreach my $id (@ids) {
    `grep -w -v $id a > a-1`;
    `cp a-1 a`;
}
`mv a-1 eval_genes.fsa`;
`cp TrainContigScores.clusters a`;
foreach my $id (@ids) {
    `grep -w -v $id a > a1`;
    `cp a1 a`;
}
`mv a1 TrainContigScores.clusters`;
`cp TrainContigScores.signals a`;
foreach my $id (@ids) {
    `grep -w -v $id a > a2`;
    `cp a2 a`;
}
`mv a2 TrainContigScores.signals`;
foreach my $type (@types) {
    foreach my $gc (@gc) {
	open(TC, "<TrainContigScores.$type.$gc") || die;
	open(OUT, ">a") || die;
	while(defined($line = <TC>)) {
	    if ($line =~ />(\S+)\:\d\:/) { 
		$id = $1; print OUT $line;
	    } 
	    else { 
		print OUT "$id $line";
	    }
	}
	close(OUT); 
	close(TC);
	foreach my $id (@ids) {
	    `grep -w -v $id a > a3`;
	    `cp a3 a`;
	}
	open(TC, "<a3") || die;
	open(OUT, ">a3.$type.$gc") || die;
	while(defined($line = <TC>)) {
	    if ($line =~ />(\S+)\:\d\:/) {
		print OUT $line;
	    } else {
		$line =~ /\S+\s+(\S+)$/; 
		print OUT "$1\n";
	    }
	}
	close(TC);
	close(OUT);
	`mv a3.$type.$gc TrainContigScores.$type.$gc`;
    }
}
`rm a a3`;
