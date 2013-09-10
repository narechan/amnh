#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:o:', \%opts);
my $file  = $opts{'f'};
my $outdir = $opts{'o'};

`mkdir -p $outdir`;

my $data = {};
my $counter = 0;
my $exoncounter = {};
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$file");
while (my $sequence_obj = $seqin->next_seq()){
    $counter++;
#    print STDERR "$counter\n";

    my $id       = $sequence_obj->display_id();
    my @id = split (/\_/, $id);
    my $seq      = $sequence_obj->seq();

    my $exon = $id[0] . "_" . $id[1] . "_" . $id[3] . "_" . $id[4];
    my $taxon = $id[2];
    $exoncounter->{$exon}++;
    
    if (exists ($data->{$exon})){
	my $exonmod = $exon . "_2";
	$data->{$exonmod}->{$taxon} = $seq;
    }
    else{
	$data->{$exon}->{$taxon} = $seq;
    }
}

foreach my $exon (keys %$data){
    open (F, ">$outdir/$exon.fa");
    
    foreach my $taxa (sort keys %{$data->{$exon}}){
	print F ">$taxa\n$data->{$exon}->{$taxa}\n";
    }
    close (F);
}
    
