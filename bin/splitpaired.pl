#!/usr/bin/perl

use strict;

my $name;
my $fasta=0; #next line is fasta
my @sequence;
my @qual;

open(ONE, ">split1");
open(TWO, ">split2");

while (<>) {
    chomp;
    s/\r//;
    my $line = $_;
    if ($line =~ /^\@(\S+)/) {
	$name = $1;
#	chop($name);
	$fasta=1;
#	print "\@$name\n";
    } elsif ($line =~ /^\+/) {
	$fasta=0;
#	print "\+$name\n";
    } else {

	my $length = length($line);
	my $half = $length/2;

	if ($fasta == 1) {
	    $sequence[0] = substr $line, 0, $half;
	    $sequence[1] = substr $line, $half;


	} elsif ($fasta == 0) {
	    $qual[0] = substr $line, 0, $half;
	    $qual[1] = substr $line, $half;

	}


	if ($fasta == 0) {
	    print ONE "\@$name"."#0/"."1\n";
	    print ONE "$sequence[0]\n";
	    print ONE "\+$name"."#0/"."1\n";
	    print ONE "$qual[0]\n";
	    print TWO "\@$name"."#0/"."2\n";
	    print TWO "$sequence[1]\n";
	    print TWO "\+$name"."#0/"."2\n";
	    print TWO "$qual[1]\n";
	}
    }
}
