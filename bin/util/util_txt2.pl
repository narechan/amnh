#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('d:', \%opts);
my $dir  = $opts{'d'};

opendir (D, "$dir");
my @files = sort (readdir (D));
shift @files;
shift @files;
closedir (D);

foreach my $file (@files){
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$dir/$file");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id       = $sequence_obj->display_id();
	my $seq      = $sequence_obj->seq();
	$seq =~s/\-//g;
	
	if ($id eq "dmel"){
	    print ">$file-$id\n$seq\n";
	}
	else {
	    next;
	}
    }
}


