#!/usr/bin/perl -w

=head1 NAME

gff3_to_fasta.pl

=head1 SYNOPSIS

gff3_to_fasta.pl

Options:

--gff3 is your gff3 annotation file
--fasta is your contig / chrom fasta file
--outdir is your output dir

Requires the bioperl libs. 

=head1 DESCRIPTION

This program simulates reads from a source fasta file
and does so using both strands.

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

#####SETUP#####
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Index::Fasta;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $fasta, $gff3, $outdir);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'g|gff3=s'        => \$gff3,
    'o|outdir=s'      => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($fasta, $gff3, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# index the fasta file
my $index = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$index->make_index($fasta);

# get the fasta name
my $fasta_name;
if ($fasta =~/\//g){
    $fasta =~m/.*\/(.*)$/;
    $fasta_name = $1;
}

else {
    $fasta_name = $fasta;
}

$fasta_name =~s/\.fasta//g;

# parse the gff3 file
open (G, "$gff3");
open (N, ">$outdir/$fasta_name.fna");
#open (P, ">$outdir/$fasta_name.faa"); 
my $counter = 0;
while (my $line = <G>){
    chomp $line;
    $counter++;

    next if ($line =~m/^\#/);
    
    my ($chrom,
	$source,
	$type,
	$start,
	$end,
	$score,
	$strand,
	$phase,
	$attr) = split (/\t/, $line);

    my @attrs = split (/\;/, $attr);
    my $attrs = {};
    foreach my $att (@attrs){
	my ($key, $value) = split (/\=/, $att);
	$attrs->{$key} = $value;
    }
    
    # bail if the type is not cds
    next unless ($type eq "CDS");
    
    my $location = Bio::Location::Simple->new(-start  => $start,
					      -end    => $end,
					      -strand => $strand);
    
    my $sequence = $index->fetch($chrom);
    my $subseq   = $sequence->subseq($location);
#    my $nucobj = Bio::Seq->new(-seq => $subseq, -alphabet => 'dna');
#    my $ptobj = $nucobj->translate(-frame => 0, -codontable_id => 11);
#    my $aa = $ptobj->seq;

    # substitute M for All in the translation due to
    # alt start codons in bacteria (codon table 11)
    # also get rid of trailing *
#    substr($aa, 0, 1, 'M');
#    chop ($aa);

#    my $name = $attrs->{'ID'};
#    $name =~s/fig\|//g;
#    $name =~s/\./\_/g;

#    print N ">$fasta_name#$name\n$subseq\n";
#    print P ">$fasta_name#$name\n$aa\n";
    print N ">$fasta_name#$counter\n$subseq\n";
}

