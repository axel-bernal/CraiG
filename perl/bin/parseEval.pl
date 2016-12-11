#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

my $usage = "$0 [-splice]\n";
my $splice = shift;
my $line;
my @results = ();
my @exon_results = ();
my $exonSection = 0;
while($line = <STDIN>, defined($line)) {
    if($line =~ /Gene\s+S\S+ity\s+([\d,\.]+)\%/) {
	push(@results, $1);
    }
    elsif($line =~ /^Transcript\s+S\S+ity\s+([\d,\.]+)\%/) {
	push(@results, $1);
    }	
    elsif($line =~ /^Exon\s+S\S+ity\s+([\d,\.]+)\%/) { 
	push(@results, $1);
    }
    elsif($line =~ /^Nucleotide\s+S\S+ity\s+([\d,\.]+)\%/) { 
	push(@results, $1);
    }
    elsif($line =~ /^Exon\s*$/) {
	$exonSection = 1;
    }
    elsif($line =~ /Correct\s+S\S+ity\s+([\d,\.]+)\%/) { 
	my $num = $1;
	if(!defined($splice)) {
	    if($exonSection > 2 && $exonSection <= 10) {
		unshift(@exon_results, $num);
	    }
	}
	elsif($splice eq "-splice") {
	    if($exonSection > 40) {
		unshift(@exon_results, $num);
	    }
	}
	$exonSection++;
    }
}

for(my $i = 0; $i<= $#exon_results; $i++) {
    $exon_results[$i] = sprintf("%.1f", $exon_results[$i]);
}
for(my $i = 0; $i<= $#results; $i++) {
    $results[$i] = sprintf("%.1f", $results[$i]);
}

print $results[-2], ' & ', $results[-1];
for(my $i = 0; $i<= $#exon_results; $i++) {
    if($i == 2 || $i == 3) {
	print  " & ", $exon_results[$i + 4];
    }
    elsif($i == 6 || $i == 7) {
	print  " & ", $exon_results[$i - 4];
    }
    else {
	print  " & ", $exon_results[$i];
    }
}

for(my $i = $#results - 2; $i >= 0; $i -= 2) {
    print ' & ', $results[$i - 1], ' & ', $results[$i];
}

print "\n";
