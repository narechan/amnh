#!/usr/bin/perl

# Purpose: This program uses standard blast output to harvest 
#          hit sequences.

#          Options:
#          -b is the blast report

# use the following modules                                                 
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Location::Simple;
use Bio::DB::GenBank;
use strict;

# get inputs                                                                    
my %opts = ();
getopts ('b:', \%opts);
my $blastreport    = $opts{'b'};

# read the parse file and extract the sequence
open (F, "$blastreport");
while (my $line = <F>){
    chomp $line;
    
    # parse the blast data
    my @line  = split (/\t/, $line);
    my $query = $line[0];
    my $hit   = $line[2];
    my $hsp   = $line[5];
    my $start = $line[14];
    my $end   = $line[15];
    my $ori   = $line[18];
    next if ($hit eq "No hits found");
        
    # get the hit sequence
    my $gb = new Bio::DB::GenBank(-format => 'Fasta');
    my $seq_obj = $gb->get_Seq_by_acc($hit);
    my $id = $seq_obj->id();
    my $desc = $seq_obj->desc();
    my $seq = $seq_obj->seq();
    
    # create a location object and extract the subseq
    my $location = Bio::Location::Simple->new(-start  => $start,
					      -end    => $end,
					      -strand => $ori);
    
    my $subseq   = $seq_obj->subseq($location);
    
    print ">$id\t$query\t$hsp\t$start-$end:$ori\n$subseq\n";
}
close (F);
