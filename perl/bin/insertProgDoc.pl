#!/usr/bin/perl -w
($] >= 5.004) || die "version is $] -- need perl 5.004 or greater";
# -*- perl -*-

$usage = "$0 source_file namespace\n";
$name = shift || die "$usage\n";
$namespace = shift || die "$usage\n";
open(NAME, "<$name");

while(defined($_ = <NAME>) && $_ !~ /^\*\sCopyright\s\(C\)\s2007/) {
    print $_;
}

print "*   Copyright (C) 2002-2007  Axel E. Bernal (abernal\@seas.upenn.edu)\n",
    "*   \n*   This program is free software; you can redistribute it and/or\n",
    "*   modify it under the terms of the GNU General Public License as\n",
    "*   published by the Free Software Foundation; either version 2 of the\n",
    "*   License, or (at your option) any later version.\n*   \n",
    "*   This program is distributed in the hope that it will be useful, but\n",
    "*   WITHOUT ANY WARRANTY; without even the implied warranty of\n",
    "*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n",
    "*   General Public License for more details.\n*   \n",
    "*   You should have received a copy of the GNU General Public License\n",
    "*   along with this program; if not, write to the Free Software\n",
    "*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA\n",
    "*   02111-1307, USA.\n*   \n",
    "*   The GNU General Public License is contained in the file COPYING.\n*   \n";

while(defined($_ = <NAME>)) {
    print $_;
}

close(NAME);



__END__;

if($name =~ /(\S+).cpp/) {
    while(defined($_ = <NAME>) && $_ !~ /^\s*\#include/) {
        print  $_;
    }
    
    while(defined($_) && $_ =~ /^\s*\#include/) {
        print $_;
        $_ = <NAME>;
    }

    print "\n/****************************************************************************\n",
    "* $name - part of the $namespace namespace, a general purpose\n",
    "*", " " x length($name), "    linear semi-markov structure prediction library\n",
    "*\n",
    "* Copyright (C) 2007  Axel E. Bernal (abernal\@seas.upenn.edu)\n",
    "*\n",
    "****************************************************************************/\n\n";
}
elsif($name =~ /(\S+).h/) {
    $class = $1;
    print "/****************************************************************************\n",
    "* $name - part of the $namespace namespace, a general purpose\n",
    "*", " " x length($name), "    linear semi-markov structure prediction library\n",
    "*\n",
    "* Copyright (C) 2007  Axel E. Bernal (abernal\@seas.upenn.edu)\n",
    "*\n",
    "****************************************************************************/\n\n";

    while(defined($_ = <NAME>) && $_ !~ /^\s*\#include/) {
        print $_;
    }

    while(defined($_) && $_ =~ /^\s*\#include/) {
        print $_;
        $_ = <NAME>;
    }

    print "\n/****************************************************************************\n",
    "* $class\n",
    "*\n",
    "* The $class class FILL\n",
    "****************************************************************************/\n\n";
    
}

while(defined($_)) {
    print $_;
    $_ = <NAME>;

}

close(NAME);


