#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 wrapping_locations < locations\n\t\tremoves from locations, all those genes that occur in wrapping location with an upstream TIS\n";

my $wlocs = shift || die $usage;
open(WLOCS, "$wlocs") || die $usage;

my %genes = % { &GenUtils::readLocsFormat(\*STDIN, 0, 'TRAIN', undef)};
my %maps = ();
my $id;
my $xid;
my $mid;
my %discordant_genes = ();

foreach $id (keys %genes) {
    $xid = $genes{$id}->contigId().$genes{$id}->{'Strand'};
    $mid = \$maps{$xid};
    if(!defined($mid)) {
	%$$mid = ();
    }
    
    my $exons = $genes{$id}->{'Exons'};
    for(my $i = 1; $i < scalar(@$exons); $i++) {
	$$mid->{'D'.$exons->[$i - 1]->[2]} = "$id.$i.".scalar(@$exons);
	$$mid->{'A'.$exons->[$i]->[1]} = $id.".$i.".scalar(@$exons);
    }
}

my %wgenes = % { &GenUtils::readLocsFormat(\*WLOCS, 0, 'TEST', undef)};
foreach $id (keys %wgenes) {
    $xid = $wgenes{$id}->contigId().$wgenes{$id}->{'Strand'};
    $mid =$maps{$xid};
    if(!defined($mid)) {
	next;
    }

    my $exons = $wgenes{$id}->{'Exons'};
    for(my $i = 1; $i < scalar(@$exons); $i++) {
	my $don_pos = 'D'.$exons->[$i - 1]->[2];
	my $acc_pos = 'A'.$exons->[$i]->[1];
	if(defined($mid->{$don_pos})||defined($mid->{$acc_pos})) {
	    my $token = $mid->{$don_pos};
	    if(!defined($token)) {
		$token = $mid->{$acc_pos};
	    }
	    $token =~ /^(\S+).(\d+)\.(\d+)$/;
	    my $gid = $1;
	    if($i > $2 && (scalar(@$exons)/$3 > 1.0)) {
#		print "$gid $i and $2 and ", scalar(@$exons)," and $3 and $token\n";
		$discordant_genes{$gid} = 1;
	    }	    
	}
    }
}

foreach $id (keys %genes) {
    if(!defined($discordant_genes{$id})) {
	$genes{$id}->dump(\*STDOUT, undef, 1);
    }
}

