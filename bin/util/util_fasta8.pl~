#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('d:', \%opts);
my $nucdir = $opts{'d'};

opendir (D, "$nucdir");
my @files = sort (readdir (D));
shift @files;
shift @files;
closedir (D);

foreach my $file (@files){
    my $totlen = 0;
    my $seqinnuc = Bio::SeqIO->new (-format=>'Fasta', -file=>"$nucdir/$file");
    while (my $sequence_obj = $seqinnuc->next_seq()){
	$totlen += $sequence_obj->length();
    }
    print "$file\t$totlen\n";
}
