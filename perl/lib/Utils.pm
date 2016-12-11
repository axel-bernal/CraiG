package Utils;
use lib "$ENV{CRAIG_HOME}/perl/lib";
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
             pickAtRandom
             getGeoAccumProbs
             smoothHistogram
             hashToFile
             fileToHash
             fileToArray
             byScore
             groupByPattern
             execProgram
             min
             max
             newSequence
	     );

use strict;

sub pickAtRandom {
    my $accProbs = shift @_;
    my $pick = rand($accProbs->[-1]);
    for(my $i = 0; $i < scalar(@$accProbs); $i++) {
        if($pick <= $accProbs->[$i]) {
            return $i;
        }
    }
    return scalar(@$accProbs) - 1;
}

sub getGeoAccumProbs {
    my($mean, $num) = @_;
    my @accProbs = ();
    my $p = 1.0/($mean + 1);
    my $ind = 0;
    push(@accProbs, 1.0 - exp((++$ind + 1.0)*log(1.0 - $p)));
    for(my $i = 1; $i < $num;  $i++) {
        push(@accProbs, 1.0 - exp((++$ind + 1.0)*log(1.0 - $p)));
    }
    return \@accProbs;
}

sub smoothHistogram {
    my $histgram = @_;
    my $tmpFile = "$ENV{CRAIG_HOME}/tmp/$$.hist";
    my @vals = sort {$a <=> $b} (keys % { $histgram});
    my $min = $vals[0];
    my $max = $vals[scalar(@vals) -1];
    my $hashRes;
    &hashToFile($histgram,"$tmpFile.in"," ");
    open(FILE, ">$tmpFile") || die "Can't open $tmpFile\n";
    print FILE "set output $tmpFile.out\n";
    print FILE "set term table\n";
    print FILE "set samples ",$max - $min + 1,"\n";
    print FILE "plot $tmpFile.in s c\n";
    close(FILE);
    `gnuplot $tmpFile`;
    unlink($tmpFile);
    unlink("$tmpFile.in");
    $hashRes = &fileToHash("$tmpFile.out"," ", 1);
    unlink("$tmpFile.out");
    return $hashRes;
}

sub fileToHash {
    my($file, $colSep, $numCols) = @_;
    my(%histgram) = ();
    my($line, $rowKey);
    open(FILE, "<$file") || die "Can't open $file\n";
    while(defined($line = <FILE>)) {
	if($line =~ /^\#/ || $line =~ /^\s+$/) {
	    next;
	}
	# trimming line
	$line =~ /^\s*(\S.+\S)\s*$/;
        my(@line) = split($colSep, $1);
        $rowKey = shift @line;
	if(scalar(@line) == 1 || (defined($numCols) && $numCols == 1)) {
	    $histgram{$rowKey} = $line[1];
	}
	else {
	    if(defined($numCols) && $numCols > 1) {
		for(my $i = 0; $i < $numCols; $i++) {
		    $histgram{$rowKey}->[$i] = $line[$i];
		}
	    }		
	    else {
		push(@ { $histgram{$rowKey}},@line);
	    }
	}
#	print " ", $rowKey, " ",$line[1], "\n";
    }
    close(FILE);
    return \%histgram;
}

sub hashToFile {
    my($histgram, $file, $colSep, $numCols) = @_;
    my($key) = "";
    open(FILE, ">$file") || die "Can't open file $file\n";
    foreach $key (keys %$histgram) {
	my $i;
	if(!defined($histgram->{$key})) {
	    next;
	}
	print FILE $key;
	if(defined($numCols) && $numCols > 0) {
	    if(ref($histgram->{$key}) ne "ARRAY" || $numCols == 1) {
		print FILE $colSep, $histgram->{$key}, "\n";
	    }
	    else {
		for($i = 0; $i < $numCols; $i++) {
		    print FILE $colSep,$histgram->{$key}->[$i];
		}
		print FILE "\n";
	    }
	}
	else {
	    if(ref($histgram->{$key}) eq "ARRAY") {
		print FILE $colSep, join($colSep, @ { $histgram->{$key}}), "\n";
	    }
	    else {
		print FILE $colSep, $histgram->{$key}, "\n";
	    }
	}
    }
    close(FILE);
}

sub fileToArray {
    my($file) = @_;
    my @rows = ();
    my $line;
    open(PAT, "<$file") || die "Can't open $file\n";
    while(defined($line = <PAT>)) {
	chomp($line);
	if($line =~ /^\>/) {
	    next;
	}
	# trimming line
	$line =~ /\s*(.+)\s*/;
	my @line = split(/\s+/,$1);
	if(scalar(@line) == 1) {
	    $rows[scalar(@rows)] = $line[0];
	}
	else {
	    push(@rows, \@line);
	}
    }
    return \@rows;
}

sub byScore {
    return ($b->[0] <=> $a->[0]);
}               

sub min {
    my($a, $b) = @_;
    if($a ne $b) {
        if($a > $b) {
            return $b;
        }
        else {
            return $a;
        }
    }
    else  {
        return $a;
    }
}

sub max {
    my($a, $b) = @_;
    if($a ne $b) {
        if($a > $b) {
            return $a;
        }
        else {
            return $b;
        }
    }
    else  {
        return $a;
    }                       
}

sub groupByPattern {
    my($file, $pattern) = @_;
    my %feats = ();
    my $cid;
    my $line;

    while(defined($line = <$file>)) {
        $line =~ /$pattern/;
        $cid = $2;

        if($cid =~ /^\s*$/) {
            next;
        }
        #new contig comming 
        if(!defined($feats{$cid})) {
            @{$feats{$cid}} = ();
        }
        push(@{$feats{$cid}}, $line);
    }

    return \%feats;
}

sub execProgram {
    my($program, $buffer, $cmdOpts, $file) = @_;
    if(!length($buffer)) {
        warn "size zero buffer for $cmdOpts\n";
        return;
    } 
    my $filename = "tmp.$$";
    open(TMP, ">$filename");
    print TMP join("", @{$buffer});
    print $file `$program $cmdOpts < $filename`;
    close(TMP);
    unlink($filename);
}

sub newSequence
{
    my ($length, $char) = @_;
    
    my @seq = (); # array of single-char strings
    for(my $i=0;$i<=$length;$i++)
    {
	push @seq, $char;
    }

    return \@seq;
}
1;
