#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-
use lib "$ENV{CRAIG_HOME}/perl/lib";
use Getopt::Long;
use File::Basename;
use strict;

my $usage = "usage : $0  EDIT_FILE < params\n";
my $edit_file = shift || die "$usage";
my @resources = ();
my @filters = ();
my @features = ();
my @params = ();
my %param_ovws = ();
my %feature_ovws = ();
my $line;
open(EDIT_F, "<$edit_file") || die "Can't open $edit_file";
my $section = "none";

while($line = <EDIT_F>) {
    if($line =~ /__RESOURCE_SECTION__/) {
	$section = "resources";
    }
    elsif($line =~ /__FILTER_SECTION__/) {
	$section = "filters";
    }
    elsif($line =~ /__FEATURE_SECTION__/) {
	$section = "features";
    }
    elsif($line =~ /__PARAM_SECTION__/) {
	$section = "params";
    }
    else {
	if($section eq "resources") {
	    push(@resources, $line);
	}
	elsif($section eq "filters") {
	    push(@filters, $line);
	}
	elsif($section eq "features") {
	    if($line =~ /^Overwrite\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/) {
		$feature_ovws{$1.$2.$3} = $4;
	    }
	    else {
		push(@features, $line);
	    }
	}
	elsif($section eq "params") {
	    if($line =~ /^Overwrite\s+(\S+)\s+(\d+)\s+(\d+)\s+([\d,\.\-]+)/) {
		$param_ovws{$1} = [$2, $3, $4];
	    }
	    else {
		push(@params, $line);
	    }
	}
    }
}

close(EDIT_F);    
my $state = 0;
while($_ = <STDIN>) {
    my $print_line = 1;
#    print "state ", $state, "\n";
    if($_ =~ /^\#/ || $_ !~ /\S/) {
	next;
    }
    
    if($state == 0 && $_ =~ /^Resource/) {
	$state = 1;
    }	
    elsif($state == 0 && $_ =~ /^Filter/) {
	$state = 4;
    }	
    elsif($state == 0 && $_ =~ /^Feature/) {
	$state = 6;
    }	
    elsif($state == 1) {
	if($_ =~ /^Resource/) {
	    ;
	}
	elsif($_ =~ /^\/\//) {
	    $state = 3;
	    $print_line = 0;
	}
	else {
	    $state = 2;
	}
    }
    elsif($state == 2 && $_ =~ /^\/\//) {
	$state = 1;
    }
    elsif($state == 3) {
	print join("", @resources);
	print "//\n";
	$state = 4;
    }
    elsif($state == 4 && $_ =~ /\/\//) {
	$state = 5;
	$print_line = 0;
    }
    elsif($state == 5) {
	print join("", @filters);
	print "//\n";
	$state = 6;
    }
    elsif($state == 6) {
	if($_ =~ /\/\//) {
	    print join("", @features);
	    print "//\n";
	    $print_line = 0;
	    $state = 7;
	}
	elsif($_ =~ /^(\S+)\s+(\S+)\s+(\S+)\s+.+$/) {
	    my $ovws = $feature_ovws{$1.$2.$3};
	    if(defined($ovws)) {
		$print_line = 0;
		print "$1\t$2\t$3\t$ovws\n";
	    }
	}
    }
    elsif($state == 7) {
	if($_ =~ /^(\S+)\s+(\d+)$/) {
	    my $num_subfs = $2;
	    my $ovws = $param_ovws{$1};
	    if(defined($param_ovws{$1})) {
		print $_;
#		print "defined $1 for OW $num_subfs\n";
		for(my $i = 0; $i < $num_subfs; $i++) {
		    $_ = <STDIN>;
		    print $_;
		    if($i == $ovws->[0]) {
			$_ =~ /(\d+)\s+(\d)/;
			my $subf_entries = $1;
			my $sparse = $2;
			if(!$subf_entries) {
			    next;
			}

			$print_line = 0;			
			for(my $j = 0; $j < $subf_entries; $j++) {
			    $_ = <STDIN>;
			    $_ =~ /^([^\s]+)/;
#			    print "OW $1 with \"$ovws->[2]\"\n";
			    if($sparse && $1 == $ovws->[1]) {
				print "$1\t$ovws->[2]\n";
			    }
			    elsif($j == $ovws->[1]) {
				print "$ovws->[2]\n";
			    }
			    else {
				print $_;
			    }
			}
		    }
		}
	    }
	}
    }

    if($print_line) {
	print $_;
    }
}

print join("", @resources), "\n//\n" if $state == 3;
print join("", @filters), "\n//\n" if $state == 5;
print join("", @filters), "\n//\n" if $state == 6;
print join("", @params);
