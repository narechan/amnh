#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:l:', \%opts);
my $file  = $opts{'f'};
my $list = $opts{'l'};

my $taxa = {};
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    $taxa->{$line} = 1;
}

open (F, "$file");
while (my $line = <F>){
    chomp $line;
    next if ($line =~m/^6/);

    if (exists ($taxa->{$line})){
	print ">$line\n";
    }
    else {
	print "$line\n";
    }
}
close (F);
