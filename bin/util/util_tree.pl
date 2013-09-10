#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('t:', \%opts);
my $tree  = $opts{'t'};

open (T, "$tree");
while (my $line = <T>){
    chomp $line;
    my @taxa = split (/\,/g, $line);
    
    foreach my $taxa (sort @taxa){
	$taxa =~s/\(//g;
	$taxa =~s/\)//g;
	
	print "$taxa\n";
    }
}
close (T);
