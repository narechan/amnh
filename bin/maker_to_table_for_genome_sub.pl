#!/usr/bin/perl -w

#####SETUP#####
use strict;
use Getopt::Long;
use Pod::Usage;
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

# parse the gff3 file
open (G, "$gff3");
while (my $line = <G>){
    chomp $line;
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

