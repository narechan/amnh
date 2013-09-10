#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('a:b:', \%opts);
my $afile = $opts{'a'};
my $bfile = $opts{'b'};

my @aids;
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$afile");
while (my $sequence_obj = $seqin->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    
    push (@aids, $id);
}

my $
my $seqin2 = Bio::SeqIO->new (-format=>'Fasta', -file=>"$bfile");
while (my $sequence_obj = $seqin2->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();

    push (@aids, $id);
}
