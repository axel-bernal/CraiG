#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

use lib "$ENV{CRAIG_HOME}/perl/lib";  
use Organism;
use Contig;
use GenUtils;
use Getopt::Long;
use strict;
use enum qw(STRAND_FWD=0 STRAND_COMP=1);

my $contigFile;
my $usage = "usage : $0 -ctg contigFile < tags > filter\n";
my $seqExon;

&GetOptions(
            "-ctg=s" => \$contigFile,
	    );

(defined($contigFile) && length($contigFile) > 0) || die "Contig File must be specified. ".$usage;

my $org = Organism->new();

open(CONTIG, "<$contigFile") || die "file $contigFile not found\n";
my %contigs = % {$org->loadContigs(\*CONTIG)};
close(CONTIG);

my ($seqname2tags, $phases) = @{&read_graph_elements(\*STDIN)};
my ($alphabet_map_ref, $olap, $step) = new_alphabet_map();
my %file_keys = ();

foreach my $cid (keys %contigs)
{
    my $c = $contigs{$cid};
    my $taglines_ref;

    if(exists $seqname2tags->{$cid})
    {
        $taglines_ref = $seqname2tags->{$cid};
    }else
    {
        my @arr = ();
        $taglines_ref = \@arr;
    }
    
#    print "refs ", scalar(@$taglines_ref), "\n";
    # initialize 
    my @seq_fasta = ();
    push(@seq_fasta, [@{&new_sequence($c->{'Len'},".")}]);
    push(@seq_fasta, [@{&new_sequence($c->{'Len'},".")}]);
    
    map_tagto_seq($olap, \@seq_fasta, $taglines_ref,$c,$alphabet_map_ref, $phases->{$cid}, $step, 0);
    map_tagto_seq($olap, \@seq_fasta, $taglines_ref,$c,$alphabet_map_ref, $phases->{$cid}, $step, 1);

    print_filter($cid, \*STDOUT, \@seq_fasta);
}

sub print_filter {
    my($seq, $file, $seq_fasta)  = @_;

    if(!defined($seq_fasta->[0]) &&
       !defined($seq_fasta->[1])) {
        next;
    }
    if(!defined($seq_fasta->[0])) {
        print $file ">$seq\n- x ", scalar(@{$seq_fasta->[1]}) - 1, "\t1\n";
    }
    else {
	print $file ">$seq\n";
        GenUtils::printStateArray($seq_fasta->[0], $file);
    }
    print $file "//\n";
    if(!defined($seq_fasta->[1])) {
        print $file "- x ", scalar(@{$seq_fasta->[0]}) - 1, "\t1\n";
    }
    else {
        GenUtils::printStateArray($seq_fasta->[1], $file);
    }
}


# forward strand
sub map_tagto_seq
{
    my ($olap, $seq_fasta, $taglines_ref,$contig,$alphabet_map_ref, $phases, $step, $strand)= @_;
    my $mletter;
    my $seq_len = $contig->{'Len'};
    my $c_seq;
    my $labelSet = 0;
#    print $contig->{'Id'}, " sql ",$seq_len," $strand ", scalar(@{$taglines_ref}), "\n";
    for( ; $labelSet < scalar(@{$taglines_ref}); $labelSet++) {
        my $beg = 0;
        my $this_taglines = $taglines_ref->[$labelSet]; 
#        print "hmm ",scalar @{$this_taglines}," ", $labelSet, " $strand\n";
        my $end = scalar(@{$this_taglines});
        my $inc = 1;
        my $phase = $phases->[$labelSet][$strand];     

        if($strand == 0) {
            $c_seq = $contig->{'Seq'};        
        }
        else {
            $c_seq = GenUtils::revComp($contig->{'Seq'});
            $beg = $end - 1;
            $end = -1;
            $inc = -1;
        }
        
        for(my $i = $beg; $i != $end; $i += $inc) {
            my $line = $this_taglines->[$i];
#            print "line ", $line,"\n";
            
            if($strand == 0 && $line=~/\s*NODE_INST\s+(\w+\_B)\s+/) {
                next;
            }
            
            if($strand == 1 && $line=~/\s*NODE_INST\s+(\w+\_F)\s+/) {
                next;
            }
            

            my @array = split(/\s+/,$line);
            my ($state,$start,$len) = ($array[1],$array[2],$array[3]);

#            print STDERR "hmm $line - $state $start $len $seq_len\n";
            if($strand == 1) {
                $start = $seq_len - $start - $len  + 2;
            }      

	    my $off_letter = ".";
            my $letter;
            my $numPhases = 1;

            if(isCoding($state)) {
                my $this_seq = substr($c_seq, $start - 1, $len);
                my $numStops = GenUtils::numPhaseOccurrences($this_seq, $phase, "taa", "tga", "tag");
                if($numStops > 0) {
                    print STDERR "$contig->{'Id'} $line with $numStops stops in phase $phase\t";
		    ($phase, $numStops) = GenUtils::minPhaseOccurrences($this_seq, "taa", "tga", "tag");
                    if(!$numStops) {
                        print STDERR "\tfixed with phase $phase\n";
                    }
                    else { print STDERR "\tnew phase $phase gives $numStops stops\n";}
                }
		$numPhases = 3;
		$off_letter = '-';
            }
	    elsif(isIntronic($state)) {
		$off_letter = '_';
	    }
            
            if(exists $alphabet_map_ref->{$state}) {
                $letter = $alphabet_map_ref->{$state};                               
            }
	    else {
                $letter = "!";
                print STDERR "Unknown state: $state\n";
                exit;
            }
	    
#            print STDERR $line, "  ", $start, " ", $len, " ", $start + (3 - $phase) % 3, " $phase $letter\n";
	    
	    my @segment = ();

	    my $p = $start;
	    for( ;$p<($start+$len);$p++) {
		push(@segment, $off_letter);
	    }
	    
	    $p = $start + ($numPhases - $phase) % $numPhases;
	    
            for( ; $p<($start+$len);$p += $step->{$letter}) {
		$segment[$p - $start] = $letter;
	    }

	    $p = $start;
	    for( ;$p<($start+$len);$p++) {
		my $nletter = $segment[$p - $start];
		my $eletter = $seq_fasta->[$strand][$p];

		if($eletter eq $nletter) {
                    $mletter = $nletter;
                }
                else {
                    my $pair = $eletter.$nletter;
                    $mletter = $olap->{$pair};
		    
		    if(!defined($mletter)) {
			$pair = $nletter.$eletter;
			$mletter = $olap->{$pair};
		    }
		    
                    if(!defined($mletter)) {
			if($pair =~ /[I,S,L,E,C,O,A]+/) {
			    if($pair =~ /[\_,i,l,n]/) {
				$mletter = 'A';
			    }
			    elsif($pair =~ /\-/) {
				$mletter = 'O';
			    } 
			    else {
				$mletter = 'C';
			    }
			}
                    }
		}

#            if($p == 102364 && $strand == 1 && $contig->{'Id'} eq "50724") {
#                print "hmm $line $mletter\n";
#            }
		
                if(!defined($mletter)) {
                    warn "overlap of $seq_fasta->[$strand][$p],$nletter has not been defined, assigning non igenic - if possible - as the value\n";
                    $mletter = $nletter ne '.' ? $nletter : $eletter;
                }
                
                $seq_fasta->[$strand][$p] = $mletter;
            }
	    
            if(isCoding($state)) {        
                $phase = ($phase + $len) % 3;
            }
        }
    }
}

# index from 1
sub new_sequence
{
    my ($length,$char) = @_;

    my @seq = ("."); # array of single-char strings
    for(my $i=1;$i<=$length;$i++)
    {
	push @seq, $char;
    }

    return \@seq;
}

sub read_graph_elements
{
    my ($file) = @_;

    my %seqname2ge; # seqname -> array of lines
    my $seqname="";
    my %phases;
    my $line;

    while(defined($line = <$file>))  {
        chomp $line;
        if($line=~/^>(\S+)\s*(\d*)/) { #format:>**** where **** is seqname
            $seqname = $1;
            
            if(!exists $seqname2ge{$seqname}) {
                @{$seqname2ge{$seqname}} = ();
                @{$phases{$seqname}} = ();
            }
        }
        elsif($line=~/^Label\sSet\s\d\s\d\s(\d)\s(\d)/) {
            my @phases = ($1, (3-$2)%3);
            push(@{$phases{$seqname}}, [@phases]);
            push @{$seqname2ge{$seqname}}, [()];
        }
        elsif($line=~/^NODE_INST/) {
            push @{$seqname2ge{$seqname}->[-1]}, $line; # push the line in
#	    print "name $seqname \"", scalar(@{$seqname2ge{$seqname}->[-1]}), "\"\n";
        }
    }

    return [\%seqname2ge, \%phases];
}

sub read_seq_fasta_length
{
    my ($filename) = @_;

    my %seqname2length;
    my $seqname="";
    my $length;

    open IN_FASTA,$filename or die $!;
    while(<IN_FASTA>)
    {

	chomp;
	if($_=~/^>(.*)/) #format: >**** where **** is seqname
	{

	    # record the length of the previous sequence
	    if($seqname ne "")
	    {
		$seqname2length{$seqname} = $length;
		#print "length of $seqname: $length\n";
	    }

	    $seqname = $1;
	    $length = 0;
	}else
	{
	    $length += length($_);
	}
    }

    close IN_FASTA;

    #the last seq
    if($seqname ne "")
    {
	$seqname2length{$seqname} = $length;
    }
    
    return \%seqname2length;

}

sub isCoding {
    my($state) = @_;
    return ($state =~ /EXON/);
}

sub isIntergenic {
    my($state) = @_;
    return ($state =~ /INTERGENIC/);
}

sub isIntronic {
    my($state) = @_;
    return ($state =~ /INTRON/);
}

sub new_alphabet_map
{
    my %step = (".", 1,
		"_", 1,
                "-", 1,
                "i", 1,
                "l", 1,
                "E", 3,
                "L", 3,
                "I", 3,
                "S", 3,
                "C", 3,
                );

    my %map = ("INTERGENIC_F",".",
               "INTRON_F","i",
               "EST_INTRON_F","i",
               "LONG_INTRON_F","l",
               "INIT_EXON_F","I",
               "EST_EXON_F","C",
               "INTERNAL_EXON_F","E",
               "LAST_EXON_F","L",
               "SINGLE_EXON_F","S",
               "INTERGENIC_B",".",
               "INTRON_B","i",
               "EST_INTRON_B","i",
               "LONG_INTRON_B","l",
               "INIT_EXON_B","I",
               "EST_EXON_B","C",
               "INTERNAL_EXON_B","E",
               "LAST_EXON_B","L",
               "SINGLE_EXON_B","S",
               );
    
    my %olap = (
	"..", ".",
	"._", "_",
	".i", "i",
	".l", "l",
	".-", "-",
	".E", "E",
	".I", "I",
	".S", "S",
	".L", "L",
	".C", "C",
	"_.", ".",
	"__", "_",
	"_i", "n",
	"_l", "n",
	"_-", "-",
	"i.", "i",
	"i_", "n",
	"ii", "i",
	"il", "n",
	"i-", "-",
	"l.", "l",
	"l_", "n",
	"li", "n",
	"ll", "l",
	"l-", "-",
	"n.", "n",
	"n_", "n",
	"ni", "n",
	"nl", "n",
	"n-", "-",
	"-.", "-",
	"-_", "-",
	"-i", "-",
	"-l", "-",
	"--", "-",
	);

    return (\%map, \%olap, \%step);                   
}

__END__

    my %map = ("INTERGENIC_F_0","_",
               "INTERGENIC_F_1","_",
               "INTERGENIC_F_2","_",
               "INTRON_F_0","i",
               "INTRON_F_1","j",
               "INTRON_F_2","k",
               "EST_INTRON_F_0","i",
               "EST_INTRON_F_1","i",
               "EST_INTRON_F_2","i",
               "LONG_INTRON_F_0","l",
               "LONG_INTRON_F_1","m",
               "LONG_INTRON_F_2","n",
               "INIT_EXON_F_0","I",
               "INIT_EXON_F_1","I",
               "INIT_EXON_F_2","I",
               "EST_EXON_F_0","E",
               "EST_EXON_F_1","E",
               "EST_EXON_F_2","E",
               "INTERNAL_EXON_F_0","E",
               "INTERNAL_EXON_F_1","F",
               "INTERNAL_EXON_F_2","G",
               "LAST_EXON_F_0","L",
               "LAST_EXON_F_1","M",
               "LAST_EXON_F_2","N",
               "SINGLE_EXON_F_0","S",
               "SINGLE_EXON_F_1","S",
               "SINGLE_EXON_F_2","S",
               "INTERGENIC_B_0","_",
               "INTERGENIC_B_1","_",
               "INTERGENIC_B_2","_",
               "INTRON_B_0","i",
               "INTRON_B_1","j",
               "INTRON_B_2","k",
               "LONG_INTRON_B_0","l",
               "LONG_INTRON_B_1","m",
               "LONG_INTRON_B_2","n",
               "INIT_EXON_B_0","I",
               "INIT_EXON_B_1","I",
               "INIT_EXON_B_2","I",
               "INTERNAL_EXON_B_0","E",
               "INTERNAL_EXON_B_1","F",
               "INTERNAL_EXON_B_2","G",
               "LAST_EXON_B_0","L",
               "LAST_EXON_B_1","M",
               "LAST_EXON_B_2","N",
               "SINGLE_EXON_B_0","S",
               "SINGLE_EXON_B_1","S",
               "SINGLE_EXON_B_2","S",
               );
    
    my %olap = ("-_", "_",
                "-i", "i",
                "-j", "j",
                "-k", "k",
                "-E", "E",
                "-F", "F",
                "-G", "G",
                "-l", "l",
                "-m", "m",
                "-n", "n",
                "-I", "I",
                "-S", "S",
                "-L", "L",
                "-M", "M",
                "-N", "N",
                
                "__", "_",
                "_i", "i",
                "_j", "j",
                "_k", "k",
                "_E", "E",
                "_F", "F",
                "_G", "G",
                "_l", "l",
                "_m", "m",
                "_n", "n",
                "_I", "I",
                "_S", "S",
                "_L", "L",
                "_M", "M",
                "_N", "N",
                
                "i_", "g",
                "ii", "i",
                "iE", "F",                
                
                "E_", "V",
                "Ei", "U",
                "EE", "E",
                
                "g_", "g",
                "gi", "g",
                "gE", "W",
                
                "U_", "W",
                "Ui", "U",
                "UE", "U",
                
                "Vi", "W",
                "V_", "V",
                "VE", "V",
                
                "Wi", "W",
                "W_", "W",
                "WE", "W",
                
                );

    return (\%map, \%olap, \%step);
}
