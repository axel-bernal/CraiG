#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Utils;
use GenUtils;
use strict;

$| = 1;
my $usage = "usage: $0 < ests_in_locs_fmt > clustered_ests\n";

# loading the information
my $contigs = undef;
my %ests = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TEST', undef)};
my @ests = sort {&Est::compareTo($a, $b);} (values %ests);


my $i = 0;
my $length;
my $lest = $ests[$i++];

scalar(@{$lest->{'Exons'}}) || die "$lest->{'Id'}\n";
my $est;

while($i < scalar(@ests)) {
    $est = $ests[$i];

    (scalar(@{$est->{'Exons'}}) && scalar(@{$lest->{'Exons'}})) ||
        die "$est->{'Id'} $lest->{'Id'}\n";     
    
    if($est->{'Exons'}->[0]->[0] eq $lest->{'Exons'}->[0]->[0]) {
        if($lest->highestCoord() < $est->lowestCoord()) {
            ;
        }
        else {
#            print "overlap $est->{'Id'} $lest->{'Id'}\n";
            if($est->{'Strand'} == $lest->{'Strand'}) {
                $est->{'gId'} = $lest->{'gId'};
                if($est->highestCoord() > $lest->highestCoord()) {
                    $lest = $est;
                }
            }
        }
    }
    $i++;
}

foreach $est (@ests) {
    undef $lest->{'Id'};
    $est->dump(\*STDOUT, undef, 1);
}

__END__;

&GenUtils::displayEstList(\%ests, \*STDOUT);


my $line;
my %ests = ();
my %regions = ();
my $offset;
my $tid;
my $est;
my $i;
my $j;
my $est_id;

while(defined($line = <STDIN>)) {
    if($line =~ /^\#/) {
        next;
    }
    my @line = split(/\t/, $line);

    $tid = undef;
    if($line[8] =~ /Target=Sequence\:([^\;,\s]+)/) {
        $tid = $1;
    }
    if(!defined($tid)) {
        next;
    }

    # the first element in @{$ests{$tid}} contains the total range
    if(!defined($ests{$tid})) {
        push(@{$regions{$line[0]}}, $tid);
        push(@{$ests{$tid}}, [$offset, $line[0], $line[1], $line[3], $line[4], $line[5], $line[6], $line[7]]);
#        print "head -> $tid, $offset, $line[3], $line[4]\n";
    }
    else {
        $ests{$tid}->[0]->[3] = Utils::min($ests{$tid}->[0]->[3], $line[3]);
        $ests{$tid}->[0]->[4] = Utils::max($ests{$tid}->[0]->[4], $line[4]);
#        print "head(modified) -> $tid, $offset, $ests{$tid}->[0]->[3], $ests{$tid}->[0]->[4]\n";
    }

    my $index = 1;    
    if(scalar(@{$ests{$tid}}) > 1) {
        if($line[6] eq "+") {
            while($index < scalar(@{$ests{$tid}}) && $line[3] > $ests{$tid}->[$index]->[2]) {
                $index++;
            }
        }
        else {
            while($index < scalar(@{$ests{$tid}}) && $line[4] < $ests{$tid}->[$index]->[3]) {
                $index++;            
            }
        }
    }
    splice(@{$ests{$tid}}, $index, 0, [$line[0], $line[1], $line[3], $line[4], $line[5], $line[6], $line[7]]);
#    print "tail -> $tid, ", scalar(@{$ests{$tid}}), " $line[3], $line[4]\n";
}


foreach my $region (keys %regions) {
    my @sorted_ids = @{$regions{$region}};
    my $maxRegion = -1;
    for($i = 0; $i <= $#sorted_ids; $i++) {
        if($ests{$sorted_ids[$i]}->[4] > $maxRegion) {
            $maxRegion = $ests{$sorted_ids[$i]}->[0]->[4];
        }
    }
    for($i = 0; $i <= $#sorted_ids - 1; $i++) {
        for($j = $i + 1; $j <= $#sorted_ids; $j++) {
            if($ests{$sorted_ids[$i]}->[0]->[3] > $ests{$sorted_ids[$j]}->[0]->[3]) {
                my $tmp = $sorted_ids[$i];
                $sorted_ids[$i] = $sorted_ids[$j];
                $sorted_ids[$j] = $tmp;
            }
        }
    }
    
    my %array = ();
    for($i = 0; $i <= $maxRegion; $i++) {
        push(@{$array{"+"}}, 0);
        push(@{$array{"-"}}, 0);
    }
    
    foreach $tid (@sorted_ids) {
        for($i = $ests{$tid}->[0]->[3]; $i <= $ests{$tid}->[0]->[4]; $i++) {
            $array{$ests{$tid}->[0]->[6]}->[$i] = 1;       
        }
    }

    foreach my $strand (keys %array) {
        my $interBeg = -1;
        my $interEnd = -1;    
        my $lastLocus = 0;
        my $array = $array{$strand};
        my $lastLandMark = 1;

        for($i = 1; $i <= $maxRegion; $i++) {
            if($array->[$i] == 0 && $interBeg == -1) {
                $interBeg = $i;
            }
            if($i > 1) {
                if(($array->[$i] == 1  || $i == $maxRegion) && $array->[$i - 1] == 0) {
                    $interEnd = $i - 1;
                }
            }
            
            if($interBeg != -1 && $interEnd != -1) {
#            print ">$interBeg $interEnd\n";
                $est_id = $sorted_ids[$lastLocus];
                $j = $lastLocus;
                for(; $j < scalar(@sorted_ids); $j++) {
                    $est =  $ests{$sorted_ids[$j]}->[0];
                    if($est->[6] ne $strand) {
                        next;
                    }
                    
                    ($est->[3] >= $lastLandMark) or die "$est->[1]\n";
                    
                    if($est->[4] <= $interBeg) {
                        for(my $k = 1; $k < scalar(@{$ests{$sorted_ids[$j]}}); $k++) {
                            my $est_frag = $ests{$sorted_ids[$j]}->[$k];

                            print $est_frag->[0], ":", $est->[0] + 1, "..", $est->[0] + 1000000, "\t", $est_frag->[1], "\tCDS\t", $est_frag->[2], "\t", $est_frag->[3], "\t", $est_frag->[4], "\t", $est_frag->[5], "\t", $est_frag->[6], "\test_id \"", $est_id, "\"; transcript_id \"", $sorted_ids[$j], "\";\n";
                        }
                    }
                    else {
                        last;
                    }
                }
                if($j != $lastLocus) {
#                    print "$region, $strand, $lastLandMark, $est_id ", $j - $lastLocus, "\n";
                    $lastLandMark = $interBeg + 1;
                    $lastLocus = $j;
                }
                $interBeg = -1;
                $interEnd = -1;
            }
        }
    }
}
