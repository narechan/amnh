#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('n:', \%opts);
my $nucfile = $opts{'n'};

my $counter = 1;
my $seqinnuc = Bio::SeqIO->new (-format=>'Fasta', -file=>"$nucfile");
while (my $sequence_obj = $seqinnuc->next_seq()){
#    my $id = $sequence_obj->id();
    my $seq = $sequence_obj->seq();

#    if ($id =~m/\_$/){
#	chop $id;
#    }
#    print ">$id\n$seq\n";
    print ">$counter\n$seq\n";
    $counter++;
}
