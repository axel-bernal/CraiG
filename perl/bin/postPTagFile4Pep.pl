#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $usage = "usage : $0 < tagFile > pepTagFile\n";

my $oldLine = <STDIN>;
my $line;
my $fwdpep = 0;

while(defined($line = <STDIN>)) {
    if($fwdpep) {
        if($line =~ /NODE\_INST\s(.+)/) {
            $line = "NODE_INST PEP".$1."\n";
        }
        $fwdpep = 0;
    }
    else {
        if($line =~ /EDGE\_INST\s+PEPSTART\_(\S)/) {
            my $strand = $1;
            if($strand eq 'B') {
                if($oldLine =~ /NODE_INST\s(.+)/) {
                    $oldLine = "NODE_INST PEP".$1."\n";
                }
            }
            else {
                $fwdpep = 1;
            }
        }
    }
    print $oldLine;
    $oldLine = $line;
}

print $oldLine;
