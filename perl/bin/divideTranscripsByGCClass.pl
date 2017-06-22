#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Gene;
use GeneModelProps;
use Organism;
use Contig;
use Getopt::Long;
use strict;

my $fastaFile = "";
my $contigFile = "";
my $GCClasses = "0";
my $createOffset = 0;
my $usage = "usage : $0 -createOffset -GC gcClasses -ctg contigs -F fastaFile\nOutputs fastaFile.gcClass, for all gcClass contained in gcClasses";

&GetOptions("F=s" => \$fastaFile,
	    "ctg=s" => \$contigFile,
	    "GC=s" => \$GCClasses,
	    "createOffset" => \$createOffset,
	    );

($fastaFile ne "") || die "score file must be specified. ".$usage;
(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

use enum qw(EXON=0 INTRON=1 INTERGENIC=4);

my $model = GeneModelProps->new();
my $org = Organism->new();
print STDERR "Loading sequences..\n";
open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
open(FASTA, "<$fastaFile") || die "file $fastaFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);
my %genes = % { $org->loadGenes(\*FASTA, 'TRAIN', 1)};

my %offset = ();
if($createOffset) {
    open(OFFSET, ">$contigFile.offset");
    my $offset = 0;
    foreach my $cid (keys %contigs) {
	print OFFSET "$cid\t$offset\n";
	$offset{$cid} = $offset;
	$offset += 100000000;
    }
    close(OFFSET);
}

open(OFFSET, "<$contigFile.offset") || die "Cannot open $contigFile.offset\n";
my $line = "";
while(defined($line = <OFFSET>)) {
    my @line = split(/\s+/, $line);
    $offset{$line[0]} = $line[1];
}
close(OFFSET);

my %gcClasses = %{$model->setGCClasses($GCClasses)};
my %locsByGCClass = % {$org->getTranscriptGCClasses('TRAIN', \%gcClasses)};

foreach my $class (keys %gcClasses) {
    if( -e "$fastaFile.$class") {
        unlink ("$fastaFile.$class");
    }
}
foreach my $id (keys %locsByGCClass) {
    my $class = $locsByGCClass{$id};
    open(TMP, ">>$fastaFile.$class");
    my $gene = $genes{$id};
    $gene->displayId(\*TMP);
    my $contigBeg;
    my $i;
    my $exon;

    for($i = 0; $i < scalar(@{$gene->{'Exons'}}) ; $i++) {
	$exon = $gene->{'Exons'}->[$i];
	$contigBeg = $contigs{$exon->[0]}->{'Beg'};
	print TMP $exon->[0]."_".($exon->[1] + $offset{$exon->[0]} - $contigBeg + 1)."_".($exon->[2] + $offset{$exon->[0]} - $contigBeg + 1);
	if($i < scalar(@ { $gene->{'Exons'}}) - 1) {
	    print TMP ";";
        }
    }
    
    print TMP "\n";
    close(TMP);
}

foreach my $class (keys %gcClasses) {
    if( -e "$fastaFile.$class") {
        next;
    }
    `touch "$fastaFile.$class"`;
}
close(CONTIG);
close(FASTA);
print STDERR "Done!\n";

