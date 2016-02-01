#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('l:i:', \%opts);
my $list  = $opts{'l'};
my $in    = $opts{'i'};

my $listdata = {};
open (F, "$list");
while (my $line = <F>){
    chomp $line;
    $listdata->{$line} = 1;
}
close (F);

open (I, "$in");
while (my $line = <I>){
    chomp $line;
    my @data = split (/\t/, $line);
    if (exists ($listdata->{$data[2]})){
	print "$line\n";
    }
}
close (I);
