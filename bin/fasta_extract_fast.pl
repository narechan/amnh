#!/usr/bin/perl

# Purpose: Obtain a subset of sequences from a fasta file

#          Options:
#          -f is the fasta file
#          -l is the list of accessions you want

# use the following modules                                                 
use Getopt::Std;
use strict;

# get inputs                                                                    
my %opts = ();
getopts ('l:f:', \%opts);
my $list   = $opts{'l'};
my $fasta  = $opts{'f'};

my $accs = {};
open (LIST, $list);
while (my $line = <LIST>){
    chomp $line;
    $accs->{$line} = 1;
}
close (LIST);

# iterate
my $signal = 0;
#my $counter = 0;
open (F, "$fasta");
while (my $line = <F>){
    chomp $line;
#    $counter++;
#    print STDERR "$counter\n";
    
    if ($signal == 1){
	print "$line\n";
	$signal = 0;
    }
#    if ($line =~m/^>(.*)\s.*/){ #parses illumina deflines
    if ($line =~m/^>(.*)/){
        my $id = $1;
	if (exists ($accs->{$id})){
	    $signal++;
	    print ">$id\n";
	}
    }
}
close (F);



