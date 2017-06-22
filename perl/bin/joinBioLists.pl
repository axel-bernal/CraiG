#!/usr/bin/perl -w

my $usage = "$0 type[locs/fasta] filename [-v]< list\n";

my $type = shift || die $usage;;
my $fileName = shift   || die $usage;
my $v = shift;

if(!defined($v)) {
    $v = 0;
}

my @ids = ();

while(defined($line = <STDIN>)) {
    chomp($line);
    if(length($line) > 0) {
        my @line = split(/\s+/, $line);
        push(@ids, $line[0]);
    }
}

open(TMP, $fileName);
my $id2;
my $line = <TMP>;
while(defined($line)) {
    if($line =~ /^>(\S+)\s*/) {
        $id = $1;
        $id2 = $id;
        if($type eq "locs") {
            if($line =~ /^>\S+\s+(\d\<|\<)?([^\;,\s]+)\_[^\_,\;]+\_[^\_,\;]+/) {
                $id2 = $2;
            }
        }
    }
    elsif($line =~ /^(\S+)\s*/) {
        $id2 = $1;
    }
    else {
        $line = <TMP>;
        next;
    }

    $found = 0;
    for($i = 0; $i < scalar(@ids); $i++) {
#        if($line =~ /$ids[$i]/) {
        if($id2 eq $ids[$i]) {
            $found = 1;
            last;
        }
    }
#    print "$found $id2 $i \n";
    $print = 0;
    if(($found && !$v) || (!$found && $v)) {
        $print = 1;
    }

    if($print) {
        print $line; 
        $line = <TMP>;
        if($type eq "fasta") {
            while(defined($line) && $line !~ />\S+/) {
#            $line =~ tr/actg/NNNN/;
                print $line;
                $line = <TMP>;
            }
        }
    }
    else {
        $line = <TMP>;
        if($type eq "fasta") {
            while(defined($line) && $line !~ />\S+/) {
                $line = <TMP>;
            }
        }
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
