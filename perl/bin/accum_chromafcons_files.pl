#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

my $usage = "$0 subcontigs in_dir out_dir\n";
my $subcontigs = shift || die $usage;
my $in_dir = shift || die $usage;
my $out_dir = shift || die $usage;
my @subdirs = `ls -I*cons $in_dir/`;
my $i;

open(TMP, ">$$.nl");
print TMP "\n";
close TMP;
my %cons_files = (); 
for($i = 0; $i <= $#subdirs; $i++) {
    chomp($subdirs[$i]);
    my @ssubdirs = `ls $in_dir/$subdirs[$i]/*cons`;
    for(my $j = 0; $j <= $#ssubdirs; $j++) {    
 	my $file = $ssubdirs[$j];
	chomp($file);
 	$file =~ /([^\.,\/]+)\.([^\/]+)$/;
 	$sp = $1;
 	$prefix = $2;
 	$cons_files{"$sp.$prefix"} = 1;
	`cat $file >> $out_dir/$sp.$prefix`;
	`cat $$.nl >> $out_dir/$sp.$prefix`;
    }
}

unlink "$$.nl";

open(SUBCONTIGS, "<$subcontigs") || die "$subcontigs\n";

my $line;
my %lengths = ();
while($line = <SUBCONTIGS>) {
    $line =~ /^>(\S+)\s+\S+\_(\d+)_(\d+)$/;
    $lengths{$1} = $3 - $2 + 1;
}

foreach my $file (keys %cons_files) {    
    my @file_ids = `grep ">" "$out_dir/$file" | cut -f1 | cut -f2 -d ">"`;
    my %file_ids;
    for($i = 0; $i <= $#file_ids; $i++) {
	chomp($file_ids[$i]);
	$file_ids{$file_ids[$i]} = 1;
    }
    
#    print "$file $#file_ids $file_ids[0] ", scalar(keys %lengths), "\n";
    open(CONS, ">>$out_dir/$file") || die "$file\n";
    foreach my $cid (keys %lengths) {
	if(!defined($file_ids{$cid})) {
	    print CONS ">$cid\t$lengths{$cid}\t.\n";
	}
    }
    close(CONS);
}

# checking sanity
foreach my $file (keys %cons_files) {
    print STDERR `checkExt2Sanity.pl < $out_dir/$file`;
}

#for i in `ls encodeOct05_all/`; do for j in `ls encodeOct05_all/$i/*cons`; do echo $j; done; done
