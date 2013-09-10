#!/usr/bin/perl -w

=head1 NAME

concatenator.pl

=head1 SYNOPSIS

  concatenator.pl -- this program will concatenated subaligned genes 
    and create a single matrix while tracking coordinates in a partitions block.

Options:

 --help        Show brief help and exit
 --indir       Contains the alignments to be concatenated
 --informat    Is the format of the input alignments
 --outdir      The place to dump the output files

=head1 DESCRIPTION

Usage examp:

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::AlignIO;

my ($help, $indir, $outdir, $informat);
GetOptions(
    'h|help'          => \$help,
    'i|indir=s'      => \$indir,
    'o|outdir=s'     => \$outdir,
    'f|informat=s'   => \$informat,	   
    ) or pod2usage;

pod2usage if $help;

for my $option ($indir, $outdir, $informat){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# concatenate the alignments and define charpars
my $alndata = {};
my $alnlens = {};
opendir (J, "$indir");
my @alns = sort (readdir (J));
shift @alns;
shift @alns;
closedir (J);

my $alnlentot = 0;
my %members;
foreach my $aln (@alns){
    print STDERR "$aln\n";
    my @alnname = split (/\./, $aln);
    my $seqin = Bio::SeqIO->new(-file   => "$indir/$aln",
				  -format => "$informat");
    
    my $alnlen = 0;
    while (my $seq = $seqin->next_seq()){
	my $id        = $seq->display_id;
	
#	my $idy        = $seq->display_id;
#	my @id = split (/\//, $idy);
#	my $id = $id[0];

	$members{$id}++; #assumes that all partitions have the same number and names of taxa!!

	my $sequence = $seq->seq;
	$alnlen = $seq->length;
	$alndata->{$alnname[0]}->{$id} = $sequence;
    }
    $alnlentot += $alnlen;
    $alnlens->{$alnname[0]} = $alnlen;
}

my $members = keys %members;
    
# sort print the concatenation and charpars
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$members NCHAR=$alnlentot;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

my $start = 1;
my $end;
foreach my $count (sort keys %$alndata){
    $end = $start - 1 + $alnlens->{$count};
    print CAT "[Partition $count length $alnlens->{$count} chars $start-$end]\n";
    print PRT "CHARSET $count=$start-$end;\n";
    foreach my $sp (sort keys %{ $alndata->{$count} }){
       print CAT "$sp\t$alndata->{$count}->{$sp}\n";                                                       
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
