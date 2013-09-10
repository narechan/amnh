#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:o:', \%opts);
my $file  = $opts{'f'};
my $outdir = $opts{'o'};

`mkdir -p $outdir`;

my $seqs = {};
my $lens = {};
$/ = "\n\n";
my $counter = 0;
open (F, "$file");
while (my $chunk = <F>){
    chomp $chunk;
    
    open (O, ">$outdir/$counter.fa");
    print O "$chunk\n";
    close (O);

    my $gene;
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$outdir/$counter.fa");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id       = $sequence_obj->display_id();
        my $seq      = $sequence_obj->seq();
	my $len      = $sequence_obj->length();

	my ($id1, $id2, $taxon, $exon, $totexons) = split (/\_/, $id);
	my $acc = $id1 . "_" . $id2 . "_" . $counter;
	
	if (exists ($seqs->{$acc})){
	    until (! exists($seqs->{$acc}));
	    $acc .=
	
	$seqs->{$acc}->{$taxon}->{$exon} = $seq;
	$lens->{$acc}->{$taxon}->{$exon} = $len;
    }
#    `rm $outdir/$counter.fa`;
}

foreach my $acc (sort keys %$seqs){
    open (T, ">$outdir/$acc.fasta");

    foreach my $tax (sort keys %{$seqs->{$acc}}){
	my $genestring;
	my $totlen = 0;

	foreach my $exon (sort {$a <=> $b} keys %{$seqs->{$acc}->{$tax}}){
	    $genestring .= $seqs->{$acc}->{$tax}->{$exon};
	    $totlen += $lens->{$acc}->{$tax}->{$exon};
	}
	
	my $matchstring = "-" x $totlen;
	unless ($genestring eq $matchstring){
	    print T ">$tax\n$genestring\n";
	}

    }
    close (T);
}


#	$gene = $idc[0] . "_" . $idc[1] . "_" . $idc[3] . "_" . $idc[4] . "_" . $counter;
#	my $taxa = $idc[2];

#	my $matchstring = "\\" . "-" . "{$len}";
	
#	unless ($seq =~m/$matchstring/){
#	    $seqs->{$taxa} = $seq;
#	}
#    }
    
#    open (T, ">$outdir/$gene.fasta");
#    foreach my $tax (keys %$seqs){
#	print T ">$tax\n$seqs->{$tax}\n";
#    }
#    close (T);
#    `rm $outdir/$counter.fa`;
#}
#close (F);


