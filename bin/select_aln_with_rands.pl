#!/usr/bin/perl -w

=head1 NAME

select_aln_with_rands.pl

=head1 SYNOPSIS

  select_aln_with_rands.pl -- this program will partition a larger alignment with charsets defined
    according to specifications and output the entire thing or a random subset of taxa

Options:

 --help        Show brief help and exit
 --infile      Is the input aln   
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
use TreeSupports;

my ($help, $infile, $outdir, $members, $number, $rand);
GetOptions(
    'h|help'          => \$help,
    'o|outdir=s'      => \$outdir,
    'i|infile=s'      => \$infile,
    'm|members=s'     => \$members,
    'n|number=s'      => \$number,
    'r|random'        => \$rand,
    ) or pod2usage;

pod2usage if $help;

for my $option ($outdir, $infile, $members, $number){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# parse the aln and partitions
# instantiate the object and load stuff we need                                   
my $supportobj = TreeSupports->new;
$supportobj->load_aln    ($infile);

# get the charsets                                                                
my $charsets = $supportobj->get_charsets;

# store alignment information                                                     
my $partitions = {};
my $lengths    = {};
my @taxa;
my $alnin = Bio::AlignIO->new(-file   => "$infile",
                              -format => "nexus");

# only one aln there                                                              
my $aln = $alnin->next_aln();
my $alnlen = $aln->length;
foreach my $seq ($aln->each_seq){
    my $id        = $seq->display_id;
    push (@taxa, $id);

    foreach my $charset (sort keys %$charsets){
#       print STDERR "Storing\t$id\t$charset\n";                                  

        my $coords = $charsets->{$charset};
        my ($start, $end) = split (/\-/, $coords);

        my $partition = $seq->subseq($start, $end);
	my $partlen   = length ($partition);

	my $insert;
        if ($rand){
            my @substring = split (//, $partition);
            my @shuffsubstr = shuffle @substring;
            $insert = join "", @shuffsubstr;
        }
        else{
            $insert = $partition;
	}

	$partitions->{$charset}->{$id} = $insert;
        $lengths->{$charset} = $partlen;
    }
}


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
    print CAT "DIMENSIONS NTAX=$members NCHAR=$alnlen;\n";
    print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
    print CAT "MATRIX\n";
    print PRT "BEGIN SETS;\n";

    my $startaln = 1;
    my $endaln;
    foreach my $count (sort keys %$partitions){
	$endaln = $startaln - 1 + $lengths->{$count};
	print CAT "[Partition $count length $lengths->{$count} chars $startaln-$endaln]\n";
	print PRT "CHARSET $count=$startaln-$endaln;\n";
#	foreach my $sp (sort keys %{ $alndata->{$count} }){
	foreach my $sp (sort keys %species){
	    print CAT "$sp\t$partitions->{$count}->{$sp}\n";                                             
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
