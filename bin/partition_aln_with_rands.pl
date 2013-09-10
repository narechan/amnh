#!/usr/bin/perl -w

=head1 NAME

partition_aln_with_rands.pl

=head1 SYNOPSIS

  partition_aln_with_rands.pl -- this program will partition a larger alignment
    according to specifications and output the entire thing or a random subset of taxa

Options:

 --help        Show brief help and exit
 --infile      Is the input aln   
 --chars       The number of chars in each "gene"
 --members     The number of taxa members
 --outdir      Is your output dir
 --number      Is the number of random matrices you want
 --rand        Set if you want to randomize the sequences in each partition (taxa seqs)

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
use Algorithm::Numerical::Shuffle qw /shuffle/;
use Parallel::ForkManager;

my ($help, $infile, $outdir, $chars, $members, $number, $rand);
GetOptions(
    'h|help'          => \$help,
    'o|outdir=s'      => \$outdir,
    'i|infile=s'      => \$infile,
    'x|chars=s'       => \$chars,
    'm|members=s'     => \$members,
    'n|number=s'      => \$number,
    'r|random'        => \$rand,
    ) or pod2usage;

pod2usage if $help;

for my $option ($outdir, $infile, $chars, $members, $number){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# parse the aln and build partitions                                                                         
my $alnin = Bio::AlignIO->new(-file   => "$infile",
			      -format => "nexus");
my $alnobj = $alnin->next_aln();
my $alnlen = $alnobj->length;

my $start = 1;
my $alndata = {};
my $alnlens = {};
my $counter = 0;
my $alnlentot = 0;
my %taxa;
for (my $start = 1; $start < $alnlen; $start = $start + $chars){
    $counter++;

    my $end;
    if (($start + $chars - 1) > $alnlen){
	$end = $alnlen;
    }
    else{
	$end = $start + $chars - 1;
    }

    my $alnlenlcl = $end - $start + 1;
    my $partname = 'part' . $counter;
    $alnlentot += $alnlenlcl;
    $alnlens->{$partname} = $alnlenlcl;

#    print STDERR "$partname\t$alnlenlcl\n";

    foreach my $seq ($alnobj->each_seq){
	my $id        = $seq->display_id;
	$taxa{$id} = 1;
	my $partition = $seq->subseq($start, $end);

	my $insert;
	if ($rand){
	    my @substring = split (//, $partition);
            my @shuffsubstr = shuffle @substring;
            $insert = join "", @shuffsubstr;
        }
	else{
	    $insert = $partition;
	}
	$alndata->{$partname}->{$id} = $insert;
    }
}

my @taxa = keys %taxa;

# loop for as many matrices as specified
for (my $n = 1; $n <= $number; $n++){
    `mkdir -p $outdir/$n`;

    # get a random unique subset of the taxa for this matrix
    my %species;
    for (my $j = 1; $j <= $members; $j++){
	my $tax = $taxa[int(rand(@taxa))];
	if (exists ($species{$tax})){
	    ($tax = $taxa[int(rand(@taxa))]) until (! exists($species{$tax}));
	}
	$species{$tax} = 1;
    }


    # sort print the concatenation and charpars
    open (CAT, ">$outdir/$n/matrix");
    open (PRT, ">$outdir/$n/partitions");
    print CAT "#NEXUS\n";
    print CAT "BEGIN DATA;\n";
    print CAT "DIMENSIONS NTAX=$members NCHAR=$alnlentot;\n";
    print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
    print CAT "MATRIX\n";
    print PRT "BEGIN SETS;\n";

    my $startaln = 1;
    my $endaln;
    foreach my $count (sort keys %$alndata){
	$endaln = $startaln - 1 + $alnlens->{$count};
	print CAT "[Partition $count length $alnlens->{$count} chars $startaln-$endaln]\n";
	print PRT "CHARSET $count=$startaln-$endaln;\n";
#	foreach my $sp (sort keys %{ $alndata->{$count} }){
	foreach my $sp (sort keys %species){
	    print CAT "$sp\t$alndata->{$count}->{$sp}\n";                                                       
	}
	print CAT "\n";
	$startaln = $endaln + 1;
    }

    print CAT ";\n";
    print CAT "END;\n\n";
    print PRT "END;\n";

    close (CAT);
    close (PRT);

    `cat $outdir/$n/matrix $outdir/$n/partitions > $outdir/$n/$n.nexus`;
}

#####SUBS#####
