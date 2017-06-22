#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use strict;
use Getopt::Long;
use Contig;
use Utils;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

$| = 1;
my $usage = "usage: $0 -ctg contigs -protein-align < alignment_file > locs_file\n";
my $contigFile = "";
my $prot = 0;

&GetOptions("-protein-align" => \$prot,
            "-ctg=s" => \$contigFile);

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
$/ = "\n>";
my $c = Contig->new('File' => \*CONTIG);
close(CONTIG);
$/ = "\n";    
my $line;
my $idstrand;
my $gene;
my %exons = ();

while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    chomp($line);

    my @line = split(/\s/, $line);
    my $j = 0;
    if($line[4] eq "") {
	next;
    }
    my $id = "$c->{'Id'}.$line[4]";
    my $strand = ($line[1] < $line[0]) ? STRAND_COMP : STRAND_FWD;
    if(!defined $line[6] || $line[6] eq "") {
	$line[6] = $line[5];
    }
    my $exon = [$c->{'Id'}, $line[0], $line[1], $line[2], $line[3], $line[5], $line[6]];

    if(! defined $exons{$id}) {
	my @tmp1 = ();
	my @tmp2 = ();
	$exons{$id} = [\@tmp1, \@tmp2];
    }
    
    my $exons = $exons{$id}->[$strand];

    if($strand == STRAND_COMP) {	    
	for($j = 0; $j < scalar(@{$exons}); $j++) {
	    my $last = 0;
	    if($exon->[1] > $exons->[$j]->[1] || 
	       ($exon->[1] == $exons->[$j]->[1] && $exon->[6] > $exons->[$j]->[6])) {
		$last = 1;
	    }
	    if($last) {
		last;
	    }
	}
   }
    else {
	for($j = 0; $j < scalar(@{$exons}); $j++) {
	    my $last = 0;
	    if($exon->[1] < $exons->[$j]->[1] ||
	       ($exon->[1] == $exons->[$j]->[1] && $exon->[6] > $exons->[$j]->[6])) {
		$last = 1;
	    }
	    if($last) {
		last;
	    }
	}
   }
	
    @{$exons} = (@{$exons}[0 .. $j - 1], $exon, 
		 @{$exons}[$j .. $#$exons]);
#    if($strand == STRAND_FWD) {
#	if(scalar(@$exons) && $exons->[0]->[1] > $exon->[1]) {
#	    unshift(@$exons, $exon);
#	}
#	else {
#	    push(@$exons, $exon);
#	}
 #   }
#    else {
#	if(scalar(@$exons) && $exons->[0]->[1] < $exon->[1]) {
#	    unshift(@$exons, $exon);
#	}
#	else {
#	    push(@$exons, $exon);
#	}
 #   }
}

foreach my $id (keys %exons) {
#    my $exons = $exons{$id}->[STRAND_FWD];
#    removeOlapExons($exons, STRAND_FWD);
#    $exons = $exons{$id}->[STRAND_COMP];
#    removeOlapExons($exons, STRAND_COMP);

}

my %genes = ();

foreach my $id (keys %exons) {
    my $exons = $exons{$id}->[STRAND_FWD];
    exons2Genes(\%genes, $exons, $id, STRAND_FWD);
#    print "FWD $id ", scalar(@$exons), "\n";
    $exons = $exons{$id}->[STRAND_COMP];
    exons2Genes(\%genes, $exons, $id, STRAND_COMP);
#    print "COMP $id ", scalar(@$exons), "\n";
}

foreach my $id (keys %genes) {
    my $exons = $genes{$id}->{'Exons'};
    my $hasStart = 0;
    my $hasStop = 0;
    my $hasAcceptor = 0;
    my $hasDonor = 0;
    my $strand = $genes{$id}->{'Strand'};

    if(!scalar(@$exons)) {
	next;
    }

    my $acceptor = $c->getSequence(2, $exons->[0]->[1], $exons->[0]->[1], 0, 0, $strand);
    if(length($acceptor) == 3 && ($acceptor =~ /^ag/i)) {# || $acceptor =~ /^ac/i)) {
        $hasAcceptor = 1;
    }

    my $donor = $c->getSequence(0, $exons->[-1]->[2], $exons->[-1]->[2], 2, 0, $strand);
    if(length($donor) == 3 && ($donor =~ /gt$/i)) {# || $donor =~ /gc$/i || $donor =~ /at$/i))  {
        $hasDonor = 1;
    }
    
    my $start = $c->getSequence(0, $exons->[0]->[1], $exons->[0]->[1], 2, 0, $strand);
    if(length($start) == 3 && $start =~ /atg/i) {
        $hasStart = 1;
    }

    my $stop = $c->getSequence(2, $exons->[-1]->[2], $exons->[-1]->[2], 3, 0, $strand);
    if(length($stop) == 6 && ($stop =~ /taa$/i || $stop =~ /tga$/i || $stop =~ /tag$/i)) {
        $hasStop = 1;
    }

    my $acceptStart = 0;
    my $acceptStop = 0;
 
    if($hasStart && $prot) {
	if(($exons->[0]->[4] > $exons->[0]->[3] && $exons->[0]->[3] == 1) || 
	   ($exons->[0]->[4] < $exons->[0]->[3] && $exons->[0]->[4] == 1)) {
	    $genes{$id}->{'NoStart'} = 0;
	}
    }
    
    if($hasStop) {
        if(!$hasDonor && $prot) {
	    $genes{$id}->{'NoStop'} = 0;
#	    print "$id $stop $exons->[-1]->[1] $exons->[-1]->[2] $strand\n";
            if($stop =~ /^taa/i || $stop =~ /^tga/i || $stop =~ /^tag/i) {
                if($strand == STRAND_FWD) {
                    $exons->[-1]->[2] -= 3;
                }
                else {
                    $exons->[-1]->[2] += 3;
                }
            }
        }
    }

    $genes{$id}->displayHeader(\*STDOUT, undef);
    print "\t";
    for(my $j = 0; $j < scalar(@{$genes{$id}->{'Exons'}}); $j++) {
	my $exon = $genes{$id}->{'Exons'}->[$j];
	print $exon->[5];
	if($j < scalar(@{$genes{$id}->{'Exons'}}) - 1) {
	    print ";";
	}
    }
    print "\t";
    for(my $j = 0; $j < scalar(@{$genes{$id}->{'Exons'}}); $j++) {
	my $exon = $genes{$id}->{'Exons'}->[$j];
	print $exon->[6];
	if($j < scalar(@{$genes{$id}->{'Exons'}}) - 1) {
	    print ";";
	}
    }
    print "\n";
}

#    for($j = 0; $j < scalar(@{$genes{$id}->{'Exons'}}); $j++) {
#	$exon = $genes{$id}->{'Exons'}->[$j];
#	print ">$exon->[0].$id.",scalar(@{$genes{$id}->{'Exons'}}) - 1, ".$j\t";
#	if($genes{$id}->{'NoStart'} || $j != 0) {
#	    print "<";
#	}
#	print "$exon->[0]_$exon->[1]_$exon->[2]";
#	
#	if($genes{$id}->{'NoStop'} || $j != scalar(@{$genes{$id}->{'Exons'}}) - 1) {
#	    print ">";
#	}
#
#	print "\t$exon->[5]\t$exon->[6]\n";
#
#   }
#}



sub exons2Genes {
    my($genes, $exons, $oid, $strand) = @_;

    if(! scalar(@$exons)) {
	return;
    }

    my $offset = 0;
    my $j = 0;

    for( ; $j < scalar(@$exons); $j++) {
	my $id = "$oid.$strand.$offset";
	my $exon = $exons->[$j];
	my $create_gene = 0;
	my $last_exon;
	
#	print "receiving $id $exon->[1]_$exon->[2]\n";
	
	if(!defined($genes->{$id})) {
	    $create_gene = 1;
	}
	else {
	    for(my $i = 0; $i <= $offset; $i++) {
		$id = "$oid.$strand.$i";
		my $qexons = $genes->{$id}->{'Exons'};

		if(!scalar(@$qexons)) {
		    next;
		}

		$last_exon = $qexons->[-1];
		my $going_up = 0;
		if(scalar(@$exons) > 1) {
		    my $sec2last_exon = $qexons->[-2];
		    if($last_exon->[3] > $sec2last_exon->[3]) {
			$going_up = 1;
		    }
		}

#		print "looking at $id $last_exon->[1]_$last_exon->[2]\n";
		if($last_exon->[3] < $last_exon->[4] && $exon->[3] > $exon->[4]) {
		    if($i == $offset) {
			$create_gene = 1;
			$offset++;
			$id = "$oid.$strand.$offset";
			last;
		    }
		    else {
			next;
		    }
		}
		elsif($last_exon->[3] > $last_exon->[4] && $exon->[3] < $exon->[4]) {
		    if($i == $offset) {		    
			$create_gene = 1;
			$offset++;
			$id = "$oid.$strand.$offset";
			last;
		    }
		    else {
			next;
		    }
		}
		else {
		    my $distance = -1;

		    if($going_up && $exon->[3] >= $last_exon->[3]) {
			$distance = Utils::min($exon->[3], $exon->[4]) - 
			    Utils::max($last_exon->[3], $last_exon->[4]);
		    }
		    elsif(!$going_up && $exon->[3] < $last_exon->[3]) {
			$distance = Utils::min($last_exon->[3], $last_exon->[4]) - 
			    Utils::max($exon->[3], $exon->[4]);
			
		    }		    

#		    print "$id $distance\n";
		    if($distance > 1 || $distance < 0) {
			if($i == $offset) {		    
			    $create_gene = 1;
			    $offset++;
			    $id = "$oid.$strand.$offset";
			    last;
			}
			else {
			    next;
			}
		    }
		    elsif($distance <= 1 && $distance >= 0) {
			last;
		    }
		}
	    }
	}

#	if($create_gene) {
#	    print "new gene $id $exon->[1]_$exon->[2]\n";
#	}
#	else {
#	    print  "attaching to $id $last_exon->[1]_$last_exon->[2];$exon->[1]_$exon->[2]\n";
#	}

	if($create_gene) {
	    $gene = Gene->new('Id' => $id, 'Seq' => "", 'Strand' => $strand,
			      'NoStart' => 1, 'NoStop' => 1, 'Phase' => 0);
	    $genes->{$id} = $gene;
	}

	push(@{$genes->{$id}->{'Exons'}}, $exon);
    }
}


sub removeOlapExons {
    my($exons, $strand) = @_;

    my @parts = ();

    if(! scalar(@$exons)) {
	return;
    }
    
    my $last_exon = $exons->[0];
    for(my $j = 1; $j < scalar(@$exons); $j++) {
	my $exon = $exons->[$j];
	
	if($strand == STRAND_FWD && $exon->[1] <= $last_exon->[2]) {	   
	    if(abs($exon->[2] - $exon->[1]) > abs($last_exon->[2] - $last_exon->[1])) {
		splice(@$exons, $j - 1, 1);            
		$last_exon = $exon;
	    }
	    else {
		splice(@$exons, $j, 1);            
	    }
	    $j--;
	}
	elsif($strand == STRAND_COMP && $exon->[1] >= $last_exon->[2]) {
	    
	}
	else {
	    $last_exon = $exon;
	}
    }
}

