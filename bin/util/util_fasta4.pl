#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('n:', \%opts);
my $nucfile = $opts{'n'};

#my $catseq;
my $seqinnuc = Bio::SeqIO->new (-format=>'Fasta', -file=>"$nucfile");
while (my $sequence_obj = $seqinnuc->next_seq()){
    my $id = $sequence_obj->id();
    my $seq = $sequence_obj->seq();
    print ">$id\n$seq\n";
#    $catseq .= $sequence_obj->seq(); 
}

#my ($nuc, $fa) = split (/\./, $nucfile);
#print  ">$nuc\n$catseq\n";
