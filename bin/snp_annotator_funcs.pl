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
    my @data = split (/\t/, $line);
    $listdata->{$data[1]} = $data[0];
}
close (F);

open (I, "$in");
while (my $line = <I>){
    chomp $line;
    my @data = split (/\t/, $line);
    if (exists ($listdata->{$data[10]})){
	print "$line\t$listdata->{$data[10]}\n";
    }
    else {
	print "$line\tNO FUNC CAT\n";
    }
}
close (I);
