#!/usr/bin/perl -w

my $fileName = shift || die;
my @ids = ();
while(defined($line = <STDIN>)) {
    chomp($line);
    if(length($line) > 0) {
        push(@ids, $line);
    }
}
open(TMP, $fileName);
my $id2;
while(defined($line = <TMP>)) {
    if($line =~ /^>(\S+)\s*/) {
        $id = $1;
        $id2 = $id;
#        if($line =~ /^>\S+\s+(\<\d|\<)?([^\_]+)\_\d+\_\d+/) {
#            $id2 = $2;
#        }
    }
    elsif($line =~ /^(\S+)\s*/) {
        $id2 = $1;
    }
    else {
        next;
    }

    $found = 0;
    for($i = 0; $i < scalar(@ids); $i++) {
        if($line =~ /$ids[$i]/) {
#        if($id2 eq $ids[$i]) {
            $found = 1;
            last;
        }
    }
#        print "$found $id2 $ids[$i]\n";
    if($found) {
#        print $line; 
#            $line = <TMP>;
#            $line =~ tr/actg/NNNN/;
#        print $line;
    }
    else {
            print $line; 
#            $line = <TMP>;
#            print $line;
    }
}

__END__
while(defined($_ = <STDIN>)) {
	$line1 = $_;	
	chomp($line1);
#	print "\">>>$line1\"\n";
#	$line2 = `grep "$line1" ~/Projects/CRAIG/DebugOutput/tmp/TIGR/no_comp_no_trun2.fsa`;
	$line2 = `grep -A 1 ">$line1" $fileName`;
#	$line2 = `grep "$line1" ~/Projects/CRAIG/Output/HUMAN1/genezilla_results/TIGR-GS/tigr.fsa`;
	print $line2;
}
__END__
while(defined($_ = <STDIN>)) {
    if($_ =~ /Bad Signal/) {
        $_ = <STDIN>;
        if($_ =~ /\-\-/) {
            $_ = <STDIN>;
        }	
    }	
    else {
        print $_;
    }
}	
