#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Statistics::Descriptive;

my %opts = ();
getopts ('f:', \%opts);
my $fasta  = $opts{'f'};

my @lengths;
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fasta");
while (my $sequence_obj = $seqin->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    my $len      = $sequence_obj->length();
    push (@lengths, $len);
#    open (F, ">$outdir/$id.fasta");
#    print F ">$id\n$seq\n";
#    close (F);
}

my $statobj = Statistics::Descriptive::Full->new();
$statobj->add_data(@lengths);

my $count = $statobj->count();
my $sum   = $statobj->sum();
my $mean  = $statobj->mean();

print "$fasta\t$count\t$sum\t$mean\n";
