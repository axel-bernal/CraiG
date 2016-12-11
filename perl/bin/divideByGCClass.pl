#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use SigUtils;
use Contig;
use GeneModelProps;
use Getopt::Long;
use strict;

my $fastaFile = "";
my $GCClasses = "0";
my $usage = "usage : $0 -GC gcClasses -F fastaFile";

&GetOptions(
            "GC=s" => \$GCClasses,
            "F=s" => \$fastaFile,
            );

($fastaFile ne "") || die "fasta file must be specified. ".$usage;
open(CONTIGS, "$fastaFile") || die "Cannot open $fastaFile\n";

my $c;
$/ = "\n>";

my $model = GeneModelProps->new();
my %gcClasses = %{$model->setGCClasses($GCClasses)};
my $gcClass;

while(($gcClass, undef) = each (%gcClasses)) {
    my $fName = "$fastaFile.$gcClass"; 
    if(-e $fName) {
        unlink($fName);
    }
    `touch "$fName"`;
}

while(defined($c = Contig->new('File' => \*CONTIGS))) {
    my $idClass = SigUtils::getGCClass($c->{'Seq'}, \%gcClasses);
    open(TMP, ">>$fastaFile.$idClass") || die "Can't open $fastaFile.$idClass\n";
    print TMP ">$c->{'Id'}\n$c->{'Seq'}\n";
    close(TMP);
}
      
close(CONTIGS);
$/ = "\n";
