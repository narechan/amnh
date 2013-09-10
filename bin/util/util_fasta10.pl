#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('d:', \%opts);
my $dir = $opts{'d'};

opendir (D, "$dir");
my @files = sort readdir (D);
shift @files;
shift @files;
closedir (D);

foreach my $file (@files){
    my $counter = 0;
    my $seqinnuc = Bio::SeqIO->new (-format=>'Fasta', -file=>"$dir/$file");
    while (my $sequence_obj = $seqinnuc->next_seq()){
	$counter++;
	my $id       = $sequence_obj->display_id();
	my $seq      = $sequence_obj->seq();
	print ">$file-contig$counter\n$seq\n";
    }
}

