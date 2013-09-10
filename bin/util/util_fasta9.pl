#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:', \%opts);
my $file = $opts{'f'};

# five prime location
my $flocation = Bio::Location::Simple->new(-start  => 1,
					  -end    => 505,
					  -strand => 1);
# three prime location
my $tlocation = Bio::Location::Simple->new(-start  => 7936,
                                          -end    => 8440,
					   -strand => 1);

my $seqinnuc = Bio::SeqIO->new (-format=>'Fasta', -file=>"$file");
while (my $sequence_obj = $seqinnuc->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    my $fp   = substr ($seq, 0, 504);
    my $tp   = substr ($seq, 7935, 8439);
    
    
    print ">krv3LTR\n$tp\n";
}

