#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

#for contigs
#perl -ne 'if($_ =~ />[^\|]+\|(\S+)\s+/) { print ">$1\n";} else { print $_;}' < TgondiiME49Genomic_ToxoDB-7.1.fasta  > me49.chr.fa

%genes = ();
while(<>) {
    if($_ =~ /^\#/) {
	next;
    }
    @line = split(/\s+/, $_);
    if($line[2] =~ /gene/i) {
	next;
    }
    $line[8] =~ /Parent\=(\S+)/;
    $id = $1;
    defined($id) || die "no Parent=Tag in $line[8]";
    if($id =~ /^[^\|]+\|(\S+)$/) {
	$id = $1;
    }

    if(!defined $genes{$id}) {
	@{$genes{$id}} = ();
	$chrom = $line[0];
#	if($chrom =~ /[^\|]+\|(\S+)/) {
#	    $chrom = $1;
#	}	
	push(@{$genes{$id}}, $chrom);
	push(@{$genes{$id}}, $line[6]);
    }
    
    @exon = ($line[3] - 1, $line[4]);
    my $i = 2;            
    
    while($i < scalar(@{$genes{$id}}) &&
	  $genes{$id}->[$i] < $exon[0]) {
	$i += 2;
    }
    splice(@{$genes{$id}}, $i, 0, @exon);
}

print "##name\tchrom\tstrand\texonStarts\texonEnds\n";
foreach $id (keys %genes) {
    $entry = $genes{$id};
    print "$id\t$entry->[0]\t$entry->[1]\t";
#    if($entry->[1] eq '+') {
	for($i = 2; $i < scalar(@$entry); $i += 2) {
	    print $entry->[$i];
	    if($i < scalar(@$entry) - 2) {print ",";};
	}
#    }
#    else {
#	for($i = scalar(@$entry) - 2; $i > 1; $i -= 2) {
#	    print $entry->[$i];
#	    if($i > 2) {print ",";};
#	}
#    }
    print "\t";
#    if($entry->[1] eq '+') {
	for($i = 3; $i < scalar(@$entry); $i += 2) {
	    print $entry->[$i];
	    if($i < scalar(@$entry) - 1) {print ",";};
	}
#    }
#    else {
#	for($i = scalar(@$entry) - 1; $i > 2; $i -= 2) {
#	    print $entry->[$i];
#	    if($i > 3) {print ",";};
#	}	
#    }
    print "\n";
}
