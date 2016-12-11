#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $idFile = "";
my $subsFile = "";
my $noiseFile = "";
my $noise = 0;

my $usage = "usage : $0 [Options] -F ids_to_use < file with graph elems\nWhere Options can be:\n
-S substitutions_file.\tTo substitute graph elements in all given ids
-N noise_file.\t\tTo specify noisy features for all given ids
-n noise.\t\tThe actual noise value for all features in noise_file\n";

&GetOptions("-S=s" => \$subsFile,
            "-N=s" => \$noiseFile,
            "-noise=f" => \$noise,
            "-F=s" => \$idFile,
            );


my %ids = ();
my $ids = \%ids;

if(!defined($idFile) || length($idFile) == 0) {
    undef $ids;
}
else {
    open(IDS, "$idFile") || die "Can't open $idFile\n";
    my $line;
    while(defined($line = <IDS>)) {
        $line =~ /(\S+)\s*/;
        $ids{$1} = 1;
    }
}

if(length($subsFile) > 0) {
    processSubstitutions($subsFile, $ids);
}
elsif(length($noiseFile) > 0) {
    processNoisyFeatures($noiseFile, $ids, $noise);
}
exit;

sub processSubstitutions {
    my($subsFile, $ids) = @_;
    my $id;
    my %subs = ();
    my $line;
    
    open(SUBS, "$subsFile") || die "Can't open $subsFile\n";
    
    while(defined($line = <SUBS>)) {
        if($line =~ /s\/(\S+)\/(\S+)\//) {
            $subs{$1} = $2;
        }
    }
    
    while(defined($line = <STDIN>)) {
        if($line =~ /^\>(\S+)/) {
            $id = $1;
            print $line;
            next;
        }

        if(defined($ids) && !defined($ids->{$id})) {
            print $line;
            next;
        }
        
        foreach my $subKey (keys %subs) {
            if($line =~ /\s+$subKey\s+/) {
                $line =~ s/$subKey/$subs{$subKey}/;
                last;
            }
        }
        print $line;
    }
}

sub processNoisyFeatures {
    my($noiseFile, $ids, $noise) = @_;
    my $id;
    my @featNames = ();
    my $line;

    open(NOISE, "$noiseFile") || die "Can't open $noiseFile\n";
    
    while(defined($line = <NOISE>)) {
        chomp($line);
        push(@featNames, $line);
    }

    while(defined($line = <STDIN>)) {
        if($line =~ /^\>([^\.\s]+)/) {
            $id = $1;
            chomp($line);
            print $line;
            if(!defined($ids) || defined($ids->{$id})) {
                foreach my $f (@featNames) {
                    print " $f $noise";
                }
            }
            print "\n";
            next;
        }
        print $line;
    }
}
