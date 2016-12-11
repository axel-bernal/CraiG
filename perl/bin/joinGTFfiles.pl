#!/usr/bin/perl -w
my $usage = "$0 [prediction1.gtf] ... > total\n";
my $base = 0;
my $file = shift @ARGV;
my @arr = ('B','H','E');
my $set = 0;

while(defined($file) && open(TMP, "<$file")) {
    my $counter = 1;
    my $line;
    while(defined($line = <TMP>)) {	
	if($line =~ /^\#/) { print $line;}
	else {
	    my @l = split(/\s+/, $line); 
	    $l[9] = "\"".$arr[$set].substr($l[9], 1, length($l[9]) - 1); 
	    $l[11] = "\"".$arr[$set].substr($l[11], 1, length($l[11]) - 1); 
	    print $arr[$set]; 
	    for(my $i = 0; $i < 8; $i++) { print $l[$i], "\t";} 
	    for(my $i = 8; $i < 12; $i++) { print $l[$i], " ";} 
	    print "\n";
	}
	$counter++;
    }
    close(TMP);
    $set++;
    $file = shift @ARGV;
    $base += 1000*1000000;
}
