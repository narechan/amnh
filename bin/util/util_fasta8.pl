#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:l:o:', \%opts);
my $file = $opts{'f'};
my $listfile = $opts{'l'};
my $outdir = $opts{'o'};

my $list = {};
open (L, "$listfile");
while (my $line = <L>){
    chomp $line;
    $list->{$line} = 1;
}
close (L);

my $seqs = {};
my $seqinnuc = Bio::SeqIO->new (-format=>'Fasta', -file=>"$file");
while (my $sequence_obj = $seqinnuc->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $desc     = $sequence_obj->desc();
    my $seq      = $sequence_obj->seq();

    $desc =~m/\/site_id_1=(\w*)\s\//;
    my $siteid = $1;
    
    if (exists ($list->{$siteid})){
	$seqs->{$siteid}->{$id} = $seq;
    }
    else {
	next;
    }
}

foreach my $site (keys %$list){

    open (F, ">$outdir/$site");
    foreach my $id (keys %{$seqs->{$site}}){
	print F ">$id\n$seqs->{$site}->{$id}\n";
    }
    close (F);
    
}
