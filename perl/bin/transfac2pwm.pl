#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  

my $usage = "usage :$0 <  transfacMatFile > pwmFile\n";
my $line = "";

while(defined($line = <STDIN>) && $line !~ /^PO/) {;}

my @PWM = (['chars','A','C','G','T']);
my $i = 1;
while(defined($line = <STDIN>) && $line !~ /^XX/) {
    chomp($line);
    my @line = split(/\s+/, $line);
    $line[0] = $i++;
    push(@PWM, \@line);
}

for(my $j =0; $j < 5; $j++) {
    for($i = 0; $i < scalar(@PWM); $i++) {
	print $PWM[$i]->[$j]," ";
    }
    print "\n";
}
