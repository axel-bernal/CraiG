#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

my $usage = "$0 command [command options, ..] output_file\n";
my $output = pop @ARGV || die $usage;
my $cmd = join(" ", @ARGV);
($cmd =~ /^\S+/) || die $usage;
    
`qsub -p 0 -e $output.log -o $output  -q all.q -cwd -b yes -V $ENV{HOME}/Projects/devcraig/bin/$cmd`;
