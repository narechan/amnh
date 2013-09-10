#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:', \%opts);
my $file  = $opts{'f'};


my $seqobj = Bio::SeqIO->new(-file   =>"$file",
			     -format =>"fastq");
while (my $seq = $seqobj->next_seq()){
    my $id       = $seq->display_id();
    my $seq      = $seq->seq();
    print yellow;
}
