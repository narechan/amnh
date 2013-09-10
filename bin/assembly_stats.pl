#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Statistics::Descriptive;

my %opts = ();
getopts ('f:', \%opts);
my $fasta  = $opts{'f'};

my $totlen = 0;
my @lengths;
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fasta");
while (my $sequence_obj = $seqin->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    my $length   = $sequence_obj->length();

    $totlen += $length;
    push (@lengths, $length);
}

my $statobj = Statistics::Descriptive::Full->new();
$statobj->add_data(@lengths);

my $count = $statobj->count();
my $mean  = $statobj->mean();
my $median = $statobj->median();
my $max = $statobj->max();
my $min = $statobj->min();

my $halftot = $totlen / 2;
my $locallen = 0;
my $n50;
foreach $len (sort {$b<=>$a} @lengths){
    $locallen += $len;
    
    if ($locallen >= $halftot){
	$n50 = $len;
	last;
    }
    else {
	next;
    }
}

print "$count\t$mean\t$median\t$max\t$min\t$n50\t$totlen\n";
