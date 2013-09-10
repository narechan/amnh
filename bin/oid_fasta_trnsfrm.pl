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
    $id =~s/\-/\_/g;
   
#    my $desc   = $sequence_obj->desc();
#    my @desc = split (/\s/, $desc);
#    $desc[0] =~s/\.\./to/;

    my @splits = split (/\|/, $id);

    $splits[1] =~s/\.//g;
#    my @splits = split (/\./, $id);
#    $splits[2] =~s/\-/\_/g;
    my $sequence = $sequence_obj -> seq();
#    $sequence =~ tr/a-z/A-Z/;
    
#    print ">$org#$id\n$sequence\n";
    print ">$org#$splits[1]\n$sequence\n";
#    print ">$org#$desc[0]\n$sequence\n";
}

