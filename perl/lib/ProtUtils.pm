package ProtUtils;
use lib "$ENV{CRAIG_HOME}/perl/lib";
require Exporter;
@ISA = ( Exporter );
@EXPORT= qw (
             %codon
             %start
             %stop
             %comp_start
             %comp_stop
             %aa_code
             $codons
             geneticCode
             toM
             displayIdAndSeq
             );

use strict;
my %codon;
my %start;
my %stop;
my %comp_start;
my %comp_stop;
my %aa_code;
my $codons;

sub geneticCode {   
    my($type) = @_;
    
#    if ($type ne "mycoplasma" && $type ne "seleno")    {
        $stop{"TGA"}  = 1;
        $comp_stop{"TCA"} = 1;
#    }                              

    $stop{"TAG"}  = 1;
    $comp_stop{"CTA"} = 1;

    $stop{"TAA"} = 1;
    $comp_stop{"TTA"} =1;  
    
    $start{"ATG"} = 1;
    $comp_start{"CAT"} = 1;
    
    if($type ne "eukarya") {
        $start{"GTG"} = 1;
        $start{"TTG"} = 1;               
        $comp_start{"CAA"} = 1;
        $comp_start{"CAC"} = 1; 
    }
    
    $codon{"AAA"} = "K";
    $codon{"AAC"} = "N";
    $codon{"AAG"} = "K";
    $codon{"AAT"} = "N";
    $codon{"ACA"} = "T";
    $codon{"ACC"} = "T";
    $codon{"ACG"} = "T";
    $codon{"ACT"} = "T";
    $codon{"AGA"} = "R";                 
    $codon{"AGC"} = "S";
    $codon{"AGG"} = "R";
    $codon{"AGT"} = "S";
    $codon{"ATA"} = "I";
    $codon{"ATC"} = "I";
    $codon{"ATG"} = "M";
    $codon{"ATT"} = "I";
    $codon{"CAA"} = "Q";
    $codon{"CAC"} = "H";
    $codon{"CAG"} = "Q";
    $codon{"CAT"} = "H";
    $codon{"CCA"} = "P";
    $codon{"CCC"} = "P";
    $codon{"CCG"} = "P";
    $codon{"CCT"} = "P";
    $codon{"CGA"} = "R";
    $codon{"CGC"} = "R";
    $codon{"CGG"} = "R";
    $codon{"CGT"} = "R";    
    $codon{"CTA"} = "L";
    $codon{"CTC"} = "L";
    $codon{"CTG"} = "L";
    $codon{"CTT"} = "L";
    $codon{"GAA"} = "E";
    $codon{"GAC"} = "D";
    $codon{"GAG"} = "E";
    $codon{"GAT"} = "D";
    $codon{"GCA"} = "A";
    $codon{"GCC"} = "A";
    $codon{"GCG"} = "A";
    $codon{"GCT"} = "A";
    $codon{"GGA"} = "G";
    $codon{"GGC"} = "G";
    $codon{"GGG"} = "G";
    $codon{"GGT"} = "G";
    $codon{"GTA"} = "V";
    $codon{"GTC"} = "V";
    $codon{"GTG"} = "V";
    $codon{"GTT"} = "V";    
    $codon{"TAA"} = "*";
    $codon{"TAC"} = "Y";
    $codon{"TAG"} = "*";
#    if($type eq "methanosarcina") {
#      $codon{"TAG"} = "K";
#    }
    $codon{"TAT"} = "Y";
    $codon{"TCA"} = "S";
    $codon{"TCC"} = "S";
    $codon{"TCG"} = "S";
    $codon{"TCT"} = "S";
    $codon{"TGA"} = "*";
#    if($type eq "mycoplasma")    {
#        $codon{"TGA"} = "W";
#    }                       
#    elsif ($type eq  "seleno") {
#        $codon{"TGA"} = "Z";
#    }
    
    $codon{"TGC"} = "C";
    $codon{"TGG"} = "W";
    $codon{"TGT"} = "C";
    $codon{"TTA"} = "L";
    $codon{"TTC"} = "F";
    $codon{"TTG"} = "L";
    $codon{"TTT"} = "F";

    $aa_code{"A"} = "Ala\tAlanine";
    $aa_code{"C"} = "Cys\tCysteine";
    $aa_code{"D"} = "Asp\tAspartic Acid";
    $aa_code{"E"} = "Glu\tGlutamic Acid";
    $aa_code{"F"} = "Phe\tPhenylalanine";
    $aa_code{"G"} = "Gly\tGlycine";
    $aa_code{"H"} = "His\tHistidine";
    $aa_code{"I"} = "Ile\tIsoleucine";
    $aa_code{"K"} = "Lys\tLysine";
    $aa_code{"L"} = "Leu\tLeucine";
    $aa_code{"M"} = "Met\tMethionine";
    $aa_code{"N"} = "Asn\tAsparagine";
    $aa_code{"P"} = "Pro\tProline";
    $aa_code{"Q"} = "Gln\tGlutamine";
    $aa_code{"R"} = "Arg\tArginine";
    $aa_code{"S"} = "Ser\tSerine";
    $aa_code{"T"} = "Thr\tThreonine";
    $aa_code{"V"} = "Val\tValine";
    $aa_code{"W"} = "Trp\tTryptophan";
    $aa_code{"Y"} = "Tyr\tTyrosine";
    
    $codons =
       [                                        
          ["*","TAA","TAG","TGA"],
          ["A","GCA","GCC","GCG","GCT"],
          ["C","TGC","TGT"],
          ["D","GAC","GAT"],
          ["E","GAA","GAG"],
          ["F","TTC","TTT"],
          ["G","GGA","GGC","GGG","GGT"],
          ["H","CAC","CAT"],
          ["I","ATA","ATC","ATT"],
          ["K","AAA","AAG"],
          ["L","CTA","CTC","CTG","CTT","TTA","TTG"],
          ["M","ATG"],
          ["N","AAC","AAT"],
          ["P","CCA","CCC","CCG","CCT"],
          ["Q","CAA","CAG"],
          ["R","AGA","AGG","CGA","CGC","CGG","CGT"],
          ["S","AGC","AGT","TCA","TCC","TCG","TCT"],
          ["T","ACA","ACC","ACG","ACT"],
          ["V","GTA","GTC","GTG","GTT"],
          ["W","TGG"],
          ["Y","TAC","TAT"]
       ];
    return; 
}                            


sub translate {
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...NOTE: Must call '&genetic_code()' first...
#-----------------------------------------------------------------------
    my( $dnaP,$incl_stop, $nostops, $mystops) = @_;
    my( $i,$j,$ln );
    my( $x,$y );
    my( $prot ) = "";
    my( $dna  ) = $$dnaP;
    
    $ln = length($dna);
    if($incl_stop) {
        $ln -= 3;
    }
    $prot = "X" x ($ln/3);
    $dna =~ tr/a-z/A-Z/;
    ((scalar (keys %codon))) || return $prot;
    for ($i=0,$j=0; ($i < ($ln-2)); $i += 3,$j++)
    {
        $x = substr($dna,$i,3);
        $y = $codon{$x};
        if(defined($y)) {
            substr($prot,$j,1) = $y;
            if($y eq "*") {
                $nostops->{$x}++;
            }
        }
    }
    
    if($incl_stop) {
        $x = substr($dna, $i, 3);
        $y = $stop{$x};
        $mystops->{$x}++;
        substr($prot,$j,1) = (defined($y) && $y == 1) ? "*" :
            defined($codon{$x}) ? $codon{$x} : "X";
        
    }
    
    return \$prot;      
}

sub toM {
    my( $dnaP, $mystarts ) = @_;
    my( $x ) = substr( $$dnaP, 0,3 );
    my($trans);

    $x =~ tr/a-z/A-Z/;
    $mystarts->{$x}++;
    if ( $start{$x} )
    {
        return "M";
    }
    elsif ($trans = $codon{$x})
    {
        return $trans;
    }
    else
    {
        return "X";
    }
}                                 


sub displayIdAndSeq {
    my( $id, $seq, $fh ) = @_;

    if (! defined($fh) )  { $fh = \*STDOUT; }

    print $fh ">$id\n";
    &displaySeq($seq, $fh);
}

sub displaySeq {
    my ( $seq, $fh ) = @_;
    my ( $i, $n, $ln );

    if (! defined($fh) )  { $fh = \*STDOUT; }
    print $fh "$$seq\n";
}                                     
1;
