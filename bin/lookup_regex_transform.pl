#!/usr/bin/perl -w

#####SETUP#####
use strict;
use warnings;
use Chart::Gnuplot;
use Getopt::Long;
use Pod::Usage;

my ($help, $infile, $lookupfile);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'l|lookup=s'      => \$lookupfile,
	   ) or pod2usage;
pod2usage if $help;

#####MAIN#####

my $replace = {};
open (L, "$lookupfile");
while (my $line = <L>){
    chomp $line;
    my ($a, $b) = split (/\t/, $line);
    $replace->{$a} = $b;
}
close (L);

open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    my @line = split (/\t/, $line);
    my $acc = shift @line;
    foreach my $lookup (keys %$replace){
#	$line =~s/$lookup/$replace->{$lookup}/g;
	$acc =~s/$lookup/$replace->{$lookup}/g;
    }
#    print "$line\n";
    my $line2 = join "\t", @line;
    print "$acc\t$line2\n";
#    $line =~s/(@{[join '|', map { quotemeta($_) } keys %$replace]})/$replace->{$1}/g;
#    print "$line\n";
}
close (I);
