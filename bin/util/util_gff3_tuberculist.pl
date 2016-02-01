#!/usr/bin/perl -w

#####SETUP#####
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Index::Fasta;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $gff3, $tuberfile);
GetOptions(
    'h|help'          => \$help,
    'g|gff3=s'        => \$gff3,
    't|tuber=s'      => \$tuberfile,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($gff3, $tuberfile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# parse the tuber file
my $tuber = {};
open (T, "$tuberfile");
while (my $line = <T>){
    chomp $line;
    next if ($line =~m/^Functional/);
    
    my @line = split (/\t/, $line);
    $tuber->{$line[1]} = $line[0];
}
close (T);

# parse the gff3 file
open (G, "$gff3");
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

    # bail if the type is not cds                                                                          
    next unless ($type eq "CDS");

    my @attrs = split (/\;/, $attr);
    foreach my $att (@attrs){
	my ($key, $value) = split (/\=/, $att);
	if ($key eq "Name"){
	    print "$value\t$start\t$end\t$tuber->{$value}\n";
	}
	else {
	    next;
	}
    }
}


