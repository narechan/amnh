#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

foreach my $file (@ARGV){
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$file");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id       = $sequence_obj->display_id();
	my $seq      = $sequence_obj->seq();

	print ">$file-$id\n$seq\n";
    }
}
