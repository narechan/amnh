#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('l:i:', \%opts);
my $infile    = $opts{'i'};
my $lookupfile    = $opts{'l'};

my $lookup = {};
open (L, "$lookupfile");
while (my $line = <L>){
    chomp $line;
    my ($a, $b, $c) = split (/\t/, $line);
    $lookup->{$b} = $a;
}

open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    my @data = split (/\t/, $line);
    if (exists($lookup->{$data[1]})){
	print "$data[1]\t$data[0]\tCimex_lectularius#" . "$lookup->{$data[1]}\n";
    }
    else {
	print STDERR "$line\n";
    }
}
close (F);
