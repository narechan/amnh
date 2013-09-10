#!/usr/bin/perl

# Purpose: This program uses blast_forked.pl parse output to harvest 
#          hit sequences.

#          Options:
#          -d is the dir with the blast results
#          -f is the fasta file to extract from and index

# use the following modules                                                 
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Index::Fasta;
use Bio::Location::Simple;
use strict;

# get inputs                                                                    
my %opts = ();
getopts ('d:f:', \%opts);
my $dir    = $opts{'d'};
my $fasta  = $opts{'f'};

# get fasta file name
my $fasta_name;
if ($fasta =~/\//g){
    $fasta =~m/.*\/(.*)$/;
    $fasta_name = $1;
}

else {
    $fasta_name = $fasta;
}

# index the fasta file
my $index = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$index->make_index($fasta);

# read the parse file and extract the sequence
opendir (D, "$dir/parse");
my @pfiles = sort (readdir(D));
shift @pfiles;
shift @pfiles;
closedir (D);

foreach my $file (@pfiles){
    open (F, "$dir/parse/$file");
    while (my $line = <F>){
	chomp $line;
	my @line  = split (/\t/, $line);
	my $query = $line[0];
	my $hit   = $line[2];
	my $start = $line[14];
	my $end   = $line[15];
	my $ori   = $line[18];
	
	next if ($hit eq "No hits found");
	next if ($hit eq "Length insufficient");
	
	my $location = Bio::Location::Simple->new(-start  => $start,
						  -end    => $end,
						  -strand => $ori);
    
    
	my $sequence = $index->fetch($hit);
	my $subseq   = $sequence->subseq($location);
	
	print ">$query-$fasta_name\n$subseq\n";
    }
    close (F);
}
