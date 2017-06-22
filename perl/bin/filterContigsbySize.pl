#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $minSize = 0;
my $maxSize = -1;
my $usage = "usage : $0 min_size max_size < contigFile > filtered_contigs\n";
$minSize = shift || die $usage;
$maxSize = shift || die $usage;

my $org = Organism->new();
my %contigs = % {$org->loadContigs(\*STDIN)};
my $counter = 0;
my $total = 0;
foreach my $id (keys %contigs) {
    $contigs{$id}->{'Seq'} =~ s/\s+//g;
    my $len = length($contigs{$id}->{'Seq'});
#    my $add = "-" x (1000000 - $len);
#    print "$id $len \n";
#}
#__END__
#    if($id =~ /\S+\+/) {
#        print ">$id\n", $contigs{$id}->{'Seq'}.$add."\n";
#    }
#    else {
#        print ">$id\n", $add.$contigs{$id}->{'Seq'}."\n";        
#    }

    if($len < $minSize || ($maxSize > 0 && $len > $maxSize)) {
        next;
    }
    $total += $len;
    $counter++;
    print ">$id\n", $contigs{$id}->{'Seq'}, "\n";
}

print STDERR "Average = ", $total/$counter, "\n";
