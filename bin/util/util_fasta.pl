#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('c:', \%opts);
my $cutoff  = $opts{'c'};

my $filecount = 0;
foreach my $file (@ARGV){
    $filecount++;
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$file");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id       = $sequence_obj->display_id();
	($id = 0) unless ($id);
	my $seq      = $sequence_obj->seq();
	my $len      = length ($seq);
	
	if ($len >= $cutoff){
#	    print ">$filecount-$id\n$seq\n";
	    print ">$id\n$seq\n";
	}
	else {
	    next;
	}
    }
}
