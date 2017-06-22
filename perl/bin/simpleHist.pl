#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";

#print `awk '{print $4}' < xaceptors`;
my $line = "";
my @results =  ();

while($line = <>) {
    chomp($line);
    my @line = split(/\s+/, $line);
    for(my $i = 0; $i < 4; $i++) {
	$results[$i]->{$line[$i+1]}++; 
    }
    $results[4]->{$line[2] - $line[1]}++;
    $results[5]->{$line[4] - $line[3]}++;
}

for (my $i = 0; $i < 6; $i++) { 
    foreach my $key (sort {$a <=> $b} keys %{$results[$i]}) {
	print "$key\t$results[$i]->{$key}\n";
    }
    print "\n";
}

