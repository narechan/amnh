#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my $data = {};
foreach my $file (@ARGV){
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$file");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id       = $sequence_obj->display_id();
	my $desc       = $sequence_obj->desc();
	my ($string, $gene) = split (/\+/, $desc);
	my ($st, $interval) = split (/\:/, $id);
	my $seq  = $sequence_obj->seq();
	$seq =~s/\-//g;
	$data->{$gene}->{$st} = $seq;
    }
}

foreach my $gene (keys %$data){
    open (F, ">$gene.fa");
    foreach my $st (sort keys %{$data->{$gene}}){
	print F ">$st\n$data->{$gene}->{$st}\n";
    }
    close (F);
}
