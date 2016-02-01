#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:', \%opts);
my $file = $opts{'f'};

my $seqinnuc = Bio::SeqIO->new (-format=>'Fasta', -file=>"$file");
my $counter = 0;
while (my $sequence_obj = $seqinnuc->next_seq()){
    $counter++;
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    print ">$file.$counter\n$seq\n";
}

