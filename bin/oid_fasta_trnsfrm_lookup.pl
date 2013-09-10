#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:s:o:', \%opts);
my $fasta  = $opts{'f'};
my $species    = $opts{'s'};
my $outdir = $opts{'o'};

# read fasta file
my $seqin = Bio::SeqIO -> new (-format => 'Fasta', -file => "$fasta");

my $counter = 0;
open (F, ">$outdir/$species.fasta");
open (L, ">$outdir/$species.lookup");
while (my $sequence_obj = $seqin -> next_seq()){
    $counter++;
    my $id       = $sequence_obj -> display_id();
    my $sequence = $sequence_obj -> seq();
    
    print F ">$species#$counter\n$sequence\n";
    print L "$counter\t$id\n";

}
close (F);
close (L);

