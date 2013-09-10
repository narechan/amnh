#!/usr/bin/perl

# Purpose: Obtain a subset of sequences from a fasta file

#          Options:
#          -f is the fasta file that we need indexed
#          -l is the list of accessions you want

# use the following modules                                                 
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Index::Fasta;
use strict;

# get inputs                                                                    
my %opts = ();
getopts ('l:f:', \%opts);
my $list   = $opts{'l'};
my $fasta  = $opts{'f'};

# index the fasta file
my $index = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$index->make_index($fasta);

# read in the desired list and generate the data
open (LIST, $list);
while (my $line = <LIST>){
    chomp $line;
    print STDERR "$line\n";

    my $sequence = $index->fetch($line);
    my $id       = $sequence->id();
    my $seq      = $sequence->seq();
    
    print ">$id\n$seq\n";
}
close (LIST);
