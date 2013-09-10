#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:', \%opts);
my $infile    = $opts{'i'};

open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    my @data = split (/\t/, $line);
    print ">$data[0]\n$data[1]\n";
}
close (F);

