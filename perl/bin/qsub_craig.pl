#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use Getopt::Long;
use strict;

my $usage = "$0 params input_file --dev --p=priority --o=output_file --fmt=output_fmt --b=best\n";

my $params = shift || die $usage;
my $best = "1";
my $output = "test";
my $input = shift || die;
my $fmt = "locs";
my $priority = "0";
my $dev = 0;

&GetOptions("-o=s" => \$output,
            "-dev" => \$dev,
            "-fmt=s" => \$fmt,
            "-p=s" => \$priority,
            "-best=s" => \$best,
            );
    
if($params =~ /(\S*)params(\d+)/) {
    $output = $1.$output.$2.".".$fmt;
}

my $command;
if($dev) {
    $command = "qsub -p $priority -e $output.log -o $output  -q all.q -cwd -b yes -V $ENV{HOME}/Projects/newcraig/bin/craigPredict --best=$best -m --format=$fmt $params $input";
}
else {
    $command = "qsub -p $priority -e $output.log -o $output  -q all.q -cwd -b yes -V $ENV{CRAIG_HOME}/bin/craigPredict --best=$best -m --format=$fmt $params $input";
}

print $command, "\n";
system($command);
