#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;
use File::Basename;

my $description = "gene annotation";
my $color = "black";
my $force = 0;
my $substitute = 0;
my $usage = "usage : annotation fasta source genome tag [--color COLOR --desc DESCRIPTION > GBrowse2Stubs --force --substitute]\n";

my $annotation = shift || die $usage;
my $fasta = shift || die $usage; 
my $source = shift || die $usage;
my $organism = shift || die $usage;
my $tag = shift || die $usage;

&GetOptions(
    "-color=s" => \$color,
    "-desc=s" => \$description,
    "-force" => \$force,
    "-substitute" => \$substitute,
    );

if($force) {
    unlink("$organism.lengths");
    unlink("$organism.lengths.gff3");
}

`annot2gbrowse.sh $annotation $fasta $source $organism $tag`;

my $buffer = "";
if($substitute && !-e "/etc/gbrowse2/$organism.conf") {
    # create conf file base    
    my ($contig, $contig_size) = split(/\s+/, `sort -nr $organism.lengths | head -1`);

    $buffer = "[GENERAL]\ndescription   = $organism Genome\n\ninitial landmark = $contig:1..".int($contig_size/50)."\nplugins       = FilterTest RestrictionAnnotator TrackDumper FastaDumper\nautocomplete  = 1\n\ndefault tracks = $source.$tag\n\nexamples = $contig:1..".int($contig_size/50)."\n\ndebug = 1\n\n\n[TRACK DEFAULTS]\nglyph         = gene\ndatabase      = $source.$tag\nheight        = 6\nbgcolor       = lightgrey\nfgcolor       = black\nfont2color    = blue\nlabel density = 25\nbump density  = 100\nlink          = AUTO\n";

    open(CONF, ">$organism.conf") || die "Can't open $organism.conf\n";
    print CONF $buffer;
    close CONF;
    `sudo mv $organism.conf /etc/gbrowse2`;
    `sudo chown root /etc/gbrowse2/$organism.conf`;
    `sudo chgrp root /etc/gbrowse2/$organism.conf`;
}

$buffer = "db_adaptor    = Bio::DB::SeqFeature::Store\ndb_args       = -adaptor DBI::SQLite\n                -dsn /var/lib/gbrowse2/databases/$organism/$source.$tag.sqlite\nsearch options = default +autocomplete\n";

if($substitute) {
    substitute_section("TRACK DEFAULTS", "$source.$tag:database", "/etc/gbrowse2/$organism.conf",$buffer);
}
else {
    print "This should be added to /etc/gbrowse2/$organism.conf\n";
    print "[$source.$tag:database]\n";
    print $buffer;
}

$buffer = "feature    = gene:$source\ndatabase   = $source.$tag\nglyph      = gene\nlabel      = 1\nstranded   = 1\nbgcolor    = $color\nheight     = 6\nkey        = $description\n";

if($substitute) {
    substitute_section("", "$source.$tag", "/etc/gbrowse2/$organism.conf", $buffer);
}
else {
    print "[$source.$tag]\n";
    print $buffer;
}


# cannot be run in multithreaded environment
sub substitute_section {
    my($next_section, $section, $file, $what) = @_;
    my $bfile = basename($file);
    `sudo cp $file $bfile`;
    `sudo chown abernal $bfile`;
    `sudo chgrp abernal $bfile`;

    my $buffer = "";
    open(FILE, "$bfile") || die "Can't open $bfile for reading\n";
    my $state = 0;
    my $line;
    while(defined($line = <FILE>)) {
	if($state == 0 && $line =~ /^\[$section\]$/) {
	    $buffer .= $line.$what."\n";
	    $state = 1;
	}
	elsif($state == 0 && $line =~ /^\[$next_section\]$/) { 
	    $buffer .= "[$section]\n".$what."\n";
	    $state = 2;
	}
	elsif($state == 1 && $line =~ /^\[[^\]]+\]$/) {
	    $state = 3;
	}

	print "$state $line";
	$buffer .= $line if $state != 1;
    }

    $buffer .= "\n[$section]\n".$what if $state == 0;
    close(FILE);
    open(FILE, ">$bfile") || die "Can't open file $bfile for writing\n";
    print FILE $buffer;
    close(FILE);
    `sudo mv $bfile $file`;
    `sudo chown root $file`;
    `sudo chgrp root $file`;
}
