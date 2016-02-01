#!/usr/bin/perl -w

# -i is the snp chip call file
# -o is the output dir
# -t is a taxa list (in case some alignments lack certain taxa
#    and we need to insert missing blocks)

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile, $outdir, $taxalist);
GetOptions(
	   'i|infile=s'      => \$infile,
	   'o|outdir=s'      => \$outdir,
           't|taxalist=s'    => \$taxalist,
	   );

`mkdir -p $outdir`;

#####MAIN#####

# read in the taxa list
my @taxalist;
open (F, "$taxalist");
while (my $line = <F>){
    chomp $line;
    $line =~s/\s+//g;
    push (@taxalist, $line);
}
close (F);

# parse the snp calls
my $alndata = {};
my %alltaxa;
open (I, "$infile");
while (my $line = <I>){
    chomp $line;

    my ($taxa, $gene, $string) = split (/\t/, $line);
    my @string = split (//, $string);
    my $snp    = $string[49]; #the 49th position contains the measured base
    my @gene = split (/\:/, $gene, 2);

    $taxa =~s/\s//g;
    $gene[0] =~s/\-/\_/g;

    $alltaxa{$taxa} = 1;
    if ($snp){
	push (@{$alndata->{$gene[0]}->{$taxa}}, $snp); #assumes there is a readout foreach gene(exon)/taxa pair
    }
    else {
	push (@{$alndata->{$gene[0]}->{$taxa}}, "?");
    }
}

# sort print the concatenation and charpars                                   
my $nchar = 0;
my $alnlens = {};
foreach my $gene (keys %$alndata){
    my $genechar = 0;
    foreach my $taxa (keys %{$alndata->{$gene}}){
	my @string = @{$alndata->{$gene}->{$taxa}};
	my $string = @string;
	$nchar += $string;
	$genechar = $string;
	$alnlens->{$gene} = $genechar;
	last;
    }
}
my @alltaxa = keys %alltaxa;
my $alltaxa = @alltaxa;
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$alltaxa NCHAR=$nchar;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ACTGN\" DATATYPE=DNA MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

my $start = 1;
my $end;
foreach my $gene (sort keys %$alndata){
#    print STDERR "$gene\n";
    $end = $start - 1 + $alnlens->{$gene};
    print CAT "[Partition $gene length $alnlens->{$gene} chars $start-$end]\n";
    print PRT "CHARSET $gene=$start-$end;\n";
    foreach my $taxa (sort keys %{$alndata->{$gene}}){
	my @string = @{$alndata->{$gene}->{$taxa}};
	my $string = join "", @string;
	print CAT "$taxa\t$string\n";
    }
    print CAT "\n";
    $start = $end + 1;
}

print CAT ";\n";
print CAT "END;\n\n";
print PRT "END;\n";
close (CAT);
close (PRT);

`cat $outdir/matrix $outdir/partitions > $outdir/nexus`;
