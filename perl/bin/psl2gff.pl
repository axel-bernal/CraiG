#!/usr/bin/perl
# blat2gff.pl 
# Author
#   don gilbert, 2007, gilbertd@indiana.edu

use strict;
use Getopt::Long;

my ($swap,$addmatch);
my $matchtype=undef;
my $matchpart="CDS"; # no other SO choice ?
my $source="BLAT";
my %matchids;
my $prefixmatch="_mid";
my $gffver=2;
my $percentscore= 1; # default yes?

my $optok= GetOptions(
  "swapquerytarget!",\$swap, 
  "source=s", \$source, 
  "version:i", \$gffver, 
  "matchtype:s", \$matchtype, 
  "prefixmatch=s", \$prefixmatch, 
  "percentscore!", \$percentscore,
    );

die "
usage: blat2gff < blat.psl > blat.gff
options: -source $source -matchtype EST_match -swapQueryTarget -nopercentscore -version [3|2]
" unless($optok);
  
if(defined $matchtype) { $addmatch=1; $matchtype="match" unless($matchtype); }
print "##gff-version $gffver\n";
my $ispsl=0;

while(<>) {
    $ispsl=1 ; ## if m/^psLayout/; # should do but want some slack here?
    next unless(/^\d/);
    chomp; 
    my @v= split "\t";
    if(@v==22) { shift(@v);} # ucsc psl starts with a 'bin' field?
    unless($ispsl and @v==21){ die "# error: doesnt look like psl format I know" ; }

  my( $matchscore, $mismatches, $rep_matches, $orient,
    $qid, $qsize, $qstart, $qend,
    $tid, $tsize, $tstart, $tend,
    $blocksizes, $qstarts, $tstarts,
      )= @v[0..2, 8..16, 18..20];
  
    $qstart++; $tstart++; # move to 1-origin
    if($swap) {
    ($tid, $tsize, $tstart, $tend,  $qid, $qsize, $qstart, $qend,)
	= ($qid, $qsize, $qstart, $qend, $tid, $tsize, $tstart, $tend,);
    ($qstarts,$tstarts)= ($tstarts,$qstarts);
    }
    
    my $matchid    = $qid .$prefixmatch. ++$matchids{$qid};
    my @blocksizes = split( /,/ , $blocksizes );
    my @qstarts    = split( /,/ , $qstarts );
    my @tstarts    = split( /,/ , $tstarts );
    my $npart      = @qstarts;

    $matchscore    = sprintf "%.2f",(100 * ( $matchscore ) / $qsize) if($percentscore); 
  ## BioPerl Tools/Blat.pm does this, which seems wrong:  $matchscore + $mismatches + $rep_matches 
  
    if($addmatch) {
    ## with - orient need to back-calculate Query location (gff target) 
	if($orient eq '-') { ($qstart,$qend)= ($qend,$qstart); } # print reversed for Target
    gffput( 0, $tid, $source, $matchtype, $tstart, $tend, $matchscore, $orient, ".",
	    $matchid, "$qid $qstart $qend");
	$matchscore="."; # dont duplicate in match_parts; score is for entire feature
    }
    for(my $p=0; $p<$npart; $p++) {
	my($qstart,$tstart,$len)= ($qstarts[$p], $tstarts[$p], $blocksizes[$p]);
	my $qend= $qstart+$len;
	my $tend= $tstart+$len;
	$qstart++; $tstart++; # move to 1-origin
	if($orient eq '-') { # reverse all exons for Target
	    $len   = $blocksizes[$npart-1-$p]; 
	    $qend  = $qstarts[$npart-1-$p]; # lower
	    $qstart= $qend + $len; # higher
	    $qend++;
	} 
    gffput( 1, $tid, $source, $matchpart, $tstart, $tend, $matchscore, $orient, ".",
	    $matchid, "$qid $qstart $qend");
    }
  
}

sub gffput {
    my($part,$r,$s,$t,$b,$e,$p,$o,$f,$id,$tg)=@_;
    if($gffver<3){
	my $IDtag="gene" ;
	print join("\t", $r, $s, $t, $b, $e, $p, $o, $f, "$IDtag \"$id\"; Target \"$tg\""),"\n";
    }else{
	my $IDtag=($part) ? "Parent" : "ID";
	print join("\t", $r, $s, $t, $b, $e, $p, $o, $f, "$IDtag=$id;Target=$tg"),"\n";
    }
}
