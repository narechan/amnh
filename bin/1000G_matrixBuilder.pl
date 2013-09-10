#!/usr/bin/perl -w

# This program accepts a 1000G tab snp list
# and creates a matrix from that data 

# -i is the input snp table
# -o is the outdir

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile, $outdir);
GetOptions(
	   'i|infile=s'  => \$infile,
	   'o|outdir=s'  => \$outdir,
	   );

`mkdir -p $outdir`;

#####MAIN#####

# readin the file and store data
my $snps = {};
my $taxa = 0;
my $alnlentot = 0;

open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    $taxa++;

    my ($acc, $seq) = 
	split (/\t/, $line);
    print STDERR "Working on $acc\n";
    
    # blow apart the sequence
    my @seq = split (//, $seq);
    my $seqlen = @seq;
    print STDERR "length = $seqlen\n";

    # plug in the snp chars
    my $pos = 0;
    foreach my $snp (@seq){
	$pos++;
	$snps->{$pos}->{$acc} = $snp;
    }
    $alnlentot = $pos;
}
close (I);

# sort print the matrix and charpars                                   
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$taxa NCHAR=$alnlentot;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

my $start = 1;
my $end;
foreach my $count (sort {$a <=> $b} keys %$snps){
    $end = $start;
    print CAT "[Partition snp$count]\n";
    print PRT "CHARSET snp$count=$start-$end;\n";
    foreach my $sp (sort keys %{ $snps->{$count} }){
	print CAT "$sp\t$snps->{$count}->{$sp}\n"; 
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
