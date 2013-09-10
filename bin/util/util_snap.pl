#!/usr/bin/perl                                                                                                  

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:', \%opts);
my $fasta  = $opts{'f'};

my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fasta");
while (my $sequence_obj = $seqin->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    $seq         =~s/[BDEFHIJKLMOPQRSUVWXYZx]/N/g;

    print "$id\t$seq\n";
}
