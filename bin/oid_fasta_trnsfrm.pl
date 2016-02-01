#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:o:', \%opts);
my $fasta  = $opts{'f'};
my $org    = $opts{'o'};

# read fasta file
my $seqin = Bio::SeqIO -> new (-format => 'Fasta', -file => "$fasta");

while (my $sequence_obj = $seqin -> next_seq()){
    my $id       = $sequence_obj -> display_id();
    my $desc   = $sequence_obj->desc();
    my $sequence = $sequence_obj -> seq();

#    $id =~s/[[:punct:]]/\_/g;
#    print ">$org#$id\n$sequence\n";
  
#    my @desc = split (/\s/, $desc);
#    $desc[0] =~s/\.\./to/;

    # genbank data
#    my @splits = split (/\|/, $id);
#    $splits[1] =~s/\./\_/g;
#    print ">$org#$splits[1]\n$sequence\n";
    
    # trinity data
#    my @splits = split (/\|/, $id);
#    $splits[0] =~s/\./\_/g;
#    my @desc = split (/\s/, $desc);
#    $desc[0] =~s/[[:punct:]]/_/g;
#    my $cat = $splits[0] . $desc[0];
#    print ">$org#$cat\n$sequence\n";

    # rast data
    my @splits = split (/\./, $id, 2);
    $splits[1] =~s/\.//g;
    print ">$org#$splits[1]\n$sequence\n"; 

    # bedbug data
#    my @splits = split (/\s/, $id);
#    $splits[0] =~s/\-/\_/g;
#    $splits[0] =~s/\./\_/g;
#    $splits[0] =~s/[[:punct:]]/\_/g;
#    print ">$org#$splits[0]\n$sequence\n";

#    print ">$org#$id\n$sequence\n";
   
}

