#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
my $oldId = "";
my $gc = 0;
my $init;
my $end;
my $l = "";
while(defined($l = <STDIN>)) {
    chomp($l);
    @line = split(/\s+/, $l); 
    if($oldId eq "" || $oldId ne $line[0]." ".$line[1]) {
	if($oldId ne "") {
	    $gc = $gc/$length;
	    my $gcClass = 100;
	    if($gc <= 43) { $gcClass = 43} elsif($gc <= 51) { $gcClass = 51} elsif($gc <= 57) {$gcClass = 57};
	    print $oldId, " $init $end $gcClass\n";
	}
	$gc = 0;
	$length = 0;
	$init = $line[2];
    }
    $end = $line[3];

    $gc += abs($line[3] - $line[2])*$line[4];
    $length += abs($line[3] - $line[2]);
    $oldId = $line[0]." ".$line[1];
}

$gc = $gc/$length;
my $gcClass = 100;
if($gc <= 43) { $gcClass = 43} elsif($gc <= 51) { $gcClass = 51} elsif($gc <= 57) {$gcClass = 57};
print $oldId, " $init $end $gcClass\n";
