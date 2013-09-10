#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Index::Fasta;

my %opts = ();
getopts ('f:l:o:', \%opts);
my $file  = $opts{'f'};
my $list  = $opts{'l'};
my $outdir = $opts{'o'};

`mkdir -p $outdir`;

# index the fasta file                                                                                     
my $index = Bio::Index::Fasta->new(-filename => $file . ".idx", -write_flag => 1);
$index->make_index($file);

my $counter = 0;
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    $counter++;
    print STDERR "$line\n";
    
    my ($a, $b) = split (/\t/, $line);
    
    my $aseq = $index->fetch($a);
    my $bseq = $index->fetch($b);

    my $aid  = $aseq->display_id();                                                             
    my $asq  = $aseq->seq();
    
    my $bid  = $bseq->display_id();
    my $bsq  = $bseq->seq();
    
    my @aid = split (/\_/, $aid);
    
    open (O, ">$outdir/$aid[0]_$aid[1]_$counter");
    print O ">$aid\n$asq\n>$bid\n$bsq\n";
    close (O);
}
close (L);
