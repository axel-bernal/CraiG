#!/usr/bin/perl
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
use lib "$ENV{CRAIG_HOME}/perl/lib";
use GenUtils;
use Getopt::Long;

$| = 1;
my $numParts = 10;
my $inpuFile = "";
my $fmt = "fasta";
$usage = "scramble.pl -n numParts -fmt <fasta|gtf|genbank> -i InputFile\n";

&GetOptions("-n=i" => \$numParts,
            "-i=s" => \$inputFile,
            "-fmt=s" => \$fmt,
            );

(defined($inputFile) && length($inputFile) > 0) || die "Input File must be specified. ".$usage;
open(INPUT, "<$inputFile") || die $usage;

my @IDS = @{&GenUtils::scramble($fmt, \*INPUT, $numParts)};

for ($part=0; $part<$numParts; $part++)
{
    print "\nPART: $part\n";
    for ($j=0; $j< scalar(@{$IDS[$part]}); $j++,$i++)
    {
	print "$IDS[$part]->[$j]\n";
    }
}


    
