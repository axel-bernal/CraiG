#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";

use lib "$ENV{CRAIG_HOME}/perl/lib";
use GenUtils;
use ProtUtils;

# -*- perl -*-

# usage:  translate [start] [stop] [code] [frame]  < dna_seq  > fasta_seq
#
# specifying a code of "seleno" gives translation of UGA/TGA
# to "Z" (selenocysteine issue).
#


if ( (@ARGV == 1) && ($ARGV[0] =~ m/-help/) )
{
    if  ( @help = `cat $help/translate` ) 
    { die @help; } 
    else 
    { die "no helpfile $help/translate\n"; }
}


$start = 0;
$code = "standard";
$frame = 1;
$stop = 0;

foreach $arg (@ARGV)
{
    if ($arg eq "start")
    {
	$start = 1;
    }
    elsif ($arg eq "stop")
    {
      $stop = 1;
    }
    elsif ($arg =~ /-*\d+/)
    {
      $frame = $arg;
    }
    else
    {
	$code = $arg;
    }
}

&ProtUtils::geneticCode($code);
my %mystarts = ();
my %mystops = ();
my %nostops = ();

while (defined($_ = <STDIN>) && ($_ !~ /^>/)) {}
while ($_ && ($_ =~ /^>(\S+.*)$/))
{
    $header = $1;
    chomp($header);
    my $seqframe = $frame;
    $seqframe = 4 - $1 if $header =~ /\S+\s+(\d)\<\S+/;
    @seq = ();
    while (defined($_ = <STDIN>) && ($_ !~ /^[>\/]/))
    {
        $_ =~ s/\s//g;
        push(@seq,$_);
    }
    $seq = join("",@seq);
    
    if    ($seqframe == 2)
    {
        substr($seq,0,1) = "";
    }
    elsif ($seqframe == 3)
    {
        substr($seq,0,2) = "";
    }
    elsif ($seqframe == -1)
    {
        $seq = &revComp($seq);
    }
    elsif ($seqframe == -2)
    {
	$seq =  &GenUtils::revComp($seq);
	substr($seq,0,1) = "";
    }
    elsif ($seqframe == -3)
    {
        $seq = &GenUtils::revComp($seq);
        substr($seq,0,2) = "";
    }

    $m = &ProtUtils::toM(\$seq, \%mystarts);
    $seq = &ProtUtils::translate(\$seq, $stop, \%nostops, \%mystops);
    if ($start)
    {
	substr($$seq,0,1) = $m;
    }
    $ln = length($$seq);
    &ProtUtils::displayIdAndSeq("$header ($ln aa)",$seq);
    
    while ($_ && ($_ !~ /^>/)) { $_ = <STDIN>; }
}

#print "starts\t", join("\t", %mystarts), "\n";
#print "stops\t", join("\t", %mystops), "\n";
#print "coding-stops\t", join("\t", %nostops), "\n";
