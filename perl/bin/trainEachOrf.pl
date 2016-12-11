#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
my $usage = "usage :$0 outDir < listIds\n";
my $outDir = shift || die $usage;
my $line = "";
my @gc = (43, 51, 57, 100);
#my @types = ("immcoding", "coding", "intron", "intergenic");
my @types = ("intergenic");
my @ids = ();
if (! -d $outDir) {
    `mkdir $outDir`;
}

while(defined($line = <STDIN>)) {
    chomp($line);
    push(@ids, $line);
}

# commented for working with orf set.. otherwise the code trains only one orf
#while(defined($line = <STDIN>)) {
#    my @ids = ("HSDNAMIA");
#    chomp($line);
#    print STDERR "Id: $line...\n";
    `rm a a1 a2 a3 a4`;
#    if($line eq "HSDNAMIA") {
#	next;
#    }
#    push(@ids, $line);
#    `grep -w $line < genes.fsa > $outDir/genes.fsa`;
#    foreach my $id (@ids) {
#	`grep -A 1 -w $id < train_contigs.fsa > a`;
#	`cat a >> a1`;
#    }
#    `mv a1 $outDir/train_contigs.fsa`;
#    `rm a`;
#    foreach my $id (@ids) {
#	`grep -w $id TrainContigScores.clusters > a`;
#	`cat a >> a1`;
#    }
#    `mv a1 "$outDir/TrainContigScores.clusters"`;
#
#    foreach my $id (@ids) {
#	`grep -w $id TrainContigScores.signals > a`;
#	`cat a >> a2`;
#    }
#    open(SIGIN, "<a2") || die;
#    open(SIGOUT, ">$outDir/TrainContigScores.signals") || die;
#    while(defined($line = <SIGIN>)) {
#	chomp($line);
#	my @line = split(/\s+/, $line); 
#	if($line[6] == 1) { 
##	    $line[5] = 1;
#	} 
#	else {
##	    $line[5] = -1;
#	}
#	print SIGOUT join(" ", @line),"\n";
#    }
#    close(SIGIN);
#    close(SIGOUT);
#    `rm a2`;
    foreach my $type (@types) {
	foreach my $gc (@gc) {
	    open(TC, "<TrainContigScores.$type.$gc") || die;
	    open(OUT, ">a4") || die;
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
		`grep -w $id a4 > a`;
		`cat a >> a3`;
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
	    `mv a3.$type.$gc "$outDir/TrainContigScores.$type.$gc"`;
	    `rm a3`;
	}
	`rm a`;
	`rm a4`;
    }
    # running the trainer on the new data
#    system("nohup /home/abernal/Applications/XXX-1.0/Source.restore/euk_gfinder_train -D $outDir/conf -P params");
#}
__END__

    my $returnVal = fork();
    if($returnVal == 0) {
#	system("nohup /home/abernal/Applications/XXX-1.0/Source.restore/euk_gfinder_train -D $outDir/conf -P params");
	exit();
    }
    else {
	while(1) {
	    my $line = `grep -w "Train" nohup.out | tail -1`;
	    if($line =~ /Train\s+\d+ 0 0 \d+\s+\d+/) {
		# killing the trainer
		my @job = split(/\s+/, `ps -a | grep euk_gfinder_tra | grep -v grep`);
		`kill -9 $job[0]`;
		last;
	    }
	    sleep 1;
	}
    }
}
waitpid($returnVal, 0);
