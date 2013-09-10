#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('d:o:', \%opts);
my $fastadir  = $opts{'d'};
my $outdir    = $opts{'o'};

opendir (D, "$fastadir");
my @fastas = grep (/^.+\..+$/, readdir(D));
closedir (D);

`mkdir -p $outdir`;

foreach my $fasta (@fastas){
    next unless ($fasta =~m/fasta/);
    
    open (O, ">$outdir/$fasta");
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fastadir/$fasta");
    while (my $seq = $seqin->next_seq()){
        my $id             = $seq->display_id;
	my ($acc, $coords) = split (/\//, $id);
        my $length         = $seq->length;
        my $sequence       = $seq->seq;

        my $matchstring = '\\' . '?' . "{$length}"; # a bit hokey!                                          
	unless ($sequence =~m/$matchstring/){
	    print O ">$id\n$sequence\n";
	}

    }
    close (O);
}
