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

my %locsByGCClass = % {$org->getGCClasses('TRAIN', \%gcClasses)};

foreach my $class (keys %gcClasses) {
    if( -e "$fastaFile.$class") {
        unlink ("$fastaFile.$class");
    }
}

sub displayExons {
    my ($file, $exons, $moffset) = @_;
    my $i;
    my $exon;
    
    for($i = 0; $i < scalar(@$exons) ; $i++) {
	$exon = $exons->[$i];
	print $file $exon->[0]."_".($exon->[1] + $moffset)."_".($exon->[2] + $moffset);
	if(scalar(@$exons) > 1 && $i < scalar(@$exons) - 1) {
	    print $file ";";
	}
    }
}

foreach my $class (keys %locsByGCClass) {
    my $id;
    open(TMP, ">$fastaFile.$class");
    foreach $id (@ {$locsByGCClass{$class}}) {
	my $gene = $genes{$id};
	my $contigBeg = $gene->{'Contig'}->{'Beg'};
	defined($offset{$gene->{'Contig'}->{'Id'}}) || die "no offset for $gene->{'Contig'}->{'Id'} ($id)\n";
	my $moffset = $contigBeg + $offset{$gene->{'Contig'}->{'Id'}} - 1;
	$gene->displayId(\*TMP);
	displayExons(\*TMP, $gene->{'5UTRExons'}, $moffset);
	if(scalar(@{$gene->{'5UTRExons'}}) > 0) {
	    print TMP "[";
	}
	displayExons(\*TMP, $gene->{'Exons'}, $moffset);
	if(scalar(@{$gene->{'3UTRExons'}}) > 0) {
	    print TMP "]";
	}
	displayExons(\*TMP, $gene->{'3UTRExons'}, $moffset);
	print TMP "\n";
    }
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

