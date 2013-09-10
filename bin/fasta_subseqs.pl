#!/usr/bin/perl

# Purpose: Obtain sequence slices using bioperl's indexing features

#          Options:
#          -f is the fasta file that we need indexed
#          -l is the list of coordinates we want to extract as a fasta
#          -a is an optional amount of sequence 5' of your start
#          -b is an optional amount of sequence 3' of your end

# use the following modules                                                 
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Index::Fasta;
use Bio::Location::Simple;
use strict;

# get inputs                                                                    
my %opts = ();
getopts ('l:f:a:b:', \%opts);
my $list   = $opts{'l'};
my $fasta  = $opts{'f'};
my $a      = $opts{'a'};
my $b      = $opts{'b'};

# index the fasta file
my $index = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$index->make_index($fasta);

# read in the desired list and generate the data
open (LIST, $list);
while (my $line = <LIST>){
    chomp $line;
    my ($chrom, $start, $end, $ori) =
	split (/\t/, $line);

    my $startminus;
    my $endplus;
    
    my $ori = 1;
    if ($ori == 1){
	$startminus = $start - $a;
	$endplus    = $end + $b;
    }
    else {
	$startminus = $start - $b;
	$endplus    = $end + $a;
    }

    my $location = Bio::Location::Simple->new(-start  => $startminus,
					      -end    => $endplus,
					      -strand => $ori);
    
    
    my $sequence = $index->fetch($chrom);
    my $subseq   = $sequence->subseq($location);
    
    print ">$chrom:$startminus:$endplus:$ori\n$subseq\n";
}
close (LIST);
