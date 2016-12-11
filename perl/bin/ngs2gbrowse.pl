#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";
use Organism;
use Contig;
use Getopt::Long;
use strict;
use File::Basename;

my $junct_desc = "Junctions";
my $cov_desc = "Coverage";
my $force = 0;
my $substitute = 0;
my $usage = "usage : coverage junction source genome method sample_tag [--junctions-desc JUNCT_DESCRIPTION --cov-desc COV_DESCRIPTION] --force --substitute > STUB\n";

my $coverage = shift || die $usage;
my $junctions = shift || die $usage;
my $source = shift || die $usage;
my $organism = shift || die $usage;
my $method = shift || die $usage;
my $tag = shift || die $usage;

&GetOptions(
    "-junctions-desc=s" => \$junct_desc,
    "-cov-desc=s" => \$cov_desc,
    "-force" => \$force,
    "-substitute" => \$substitute,
    );

if($force) {
    unlink("$source.$method.$tag.cov_plus.bw");
    unlink("$source.$method.$tag.cov_minus.bw");
    unlink("$source.$method.$tag.junctions.gff3");
}

if($substitute && !-e "/etc/gbrowse2/$organism.conf") {
    # create conf file base
    warn "Cannot apply --substitute option as file /etc/gbrowse2/$organism.conf does not exist\n";
    $substitute = 0;
}

`ngs2gbrowse.sh $coverage $junctions $source $organism $method.$tag`;

my $buffer = "display_name = \"Forward\"\nmethod       = $method"."_$tag\nstrand       = FWD\nsource       = $source\n";

if($substitute) {
    if(!-e "/var/lib/gbrowse2/databases/$organism/metadata.txt") {
	`sudo touch /var/lib/gbrowse2/databases/$organism/metadata.txt`;
	`sudo chown apache /var/lib/gbrowse2/databases/$organism/metadata.txt`;
	`sudo chgrp apache /var/lib/gbrowse2/databases/$organism/metadata.txt`;
    }

    substitute_section("", "$source.$method.$tag.cov_plus.bw", "/var/lib/gbrowse2/databases/$organism/metadata.txt", $buffer);
}
else {
    print "This should be added to /var/lib/gbrowse2/databases/$organism/metadata.txt\n";
    print "[$source.$method.$tag.cov_plus.bw]\n";
    print $buffer;
}

if(-e "/var/lib/gbrowse2/databases/$organism/$source.$method.$tag.cov_minus.bw") {
    $buffer = "display_name = \"Reverse\"\nmethod       = $method"."_$tag\nstrand       = COMP\nsource       = $source\n";
    if($substitute) {
	substitute_section("", "$source.$method.$tag.cov_minus.bw", "/var/lib/gbrowse2/databases/$organism/metadata.txt", $buffer);
    }
    else {
	print "[$source.$method.$tag.cov_minus.bw]\n";
	print  $buffer;
    }
}

$buffer = "db_adaptor    = Bio::DB::SeqFeature::Store\ndb_args       = -adaptor DBI::SQLite\n                -dsn /var/lib/gbrowse2/databases/$organism/$source.$method.$tag.junctions.sqlite\nsearch options = default +autocomplete\n\n";
if($substitute) {
    substitute_section("TRACK DEFAULTS", "$source.$method.$tag"."_junctions:database", "/etc/gbrowse2/$organism.conf", $buffer);
}
else {
    print "This should be added to /etc/gbrowse2/$organism.conf\n";
    print "[$source.$method.$tag"."_junctions:database]\n";
    print $buffer;
}

$buffer = "db_adaptor    = Bio::DB::BigWigSet\ndb_args       = -dir /var/lib/gbrowse2/databases/$organism\n                -feature_type summary\n";

if($substitute) {
    substitute_section("TRACK DEFAULTS", "coverage:database", "/etc/gbrowse2/$organism.conf", $buffer);
}
else {
    print "[coverage:database]\n";
    print $buffer;
}

$buffer = "feature       = remark:$source\ndatabase      = $source.$method.$tag"."_junctions\nglyph         = anchored_arrow\nbump density  = 1000 \nfgcolor       =  sub  {\n                        my \$feature = shift;\n                        if (\$feature->score > 500) {\n                           return 'brown';\n                        }\n                        elsif(\$feature->score > 100) {\n                           return 'red';\n                        }\n                        elsif(\$feature->score > 10) {\n                           return 'orange';\n                        } else {\n                           return 'yellow';\n                        }\n                      }\nlabel         = 0\nlinewidth     = sub {\n                      my \$feature = shift;\n                      if (\$feature->score > 100) {\n                         return 2;\n                      } else {\n                         return 1;\n                      }\n                    }\nkey           = $junct_desc\n\n\n";

if($substitute) {
    substitute_section("", "$source.$method.$tag"."_junctions","/etc/gbrowse2/$organism.conf", $buffer);
}
else {
    print "[$source.$method.$tag"."_junctions]\n";
    print $buffer;
}

$buffer = "database     = coverage\nfeature      = $method"."_$tag:$source\nglyph        = wiggle_xyplot\nglyph_type   = boxes\nbump_density = 1\nmax_score    = 8\nclip         = 1\nscale        = both\n";

if(-e "/var/lib/gbrowse2/databases/$organism/$source.$method.$tag.cov_minus.bw") {
    $buffer .=  "min_score    = -8\nheight       = 80\n";
}
else {
    $buffer .= "min_score    = 0\nheight       = 40\n";
}

$buffer .= "label        = 0\nbgcolor      = blue\npart_color   = sub { GBrowse::Display::colorByRnaSeq(shift) }\nbump = overlap\ncolor_series  = 1\ncolor_cycle  = red:0.8 blue:0.8\nkey          = $cov_desc\ncitation     = Log2 Coverage ($cov_desc)\n";

if($substitute) {
    substitute_section("", "$source.$method.$tag"."_coverage","/etc/gbrowse2/$organism.conf", $buffer);
}
else {
    print "[$source.$method.$tag"."_coverage]\n";
    print $buffer;
}

# cannot be run in multithreaded environment
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
