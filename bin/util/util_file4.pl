#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:l:', \%opts);
my $infile  = $opts{'i'};
my $listfile = $opts{'l'};

my $list = {};
open (L, "$listfile");
while (my $line = <L>){
    chomp $line;
    $list->{$line} = 1;
}
close (L);

open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    my ($acc, $sequence) = split (/\t/, $line);
    if (exists ($list->{$acc})){
	print "$line\n";
    }
    else {
	next;
    }
}
close (F);

    
