#!/usr/bin/perl -w

# This program accepts a 1000G tab snp list
# and screens out invariant characters

# -i is the input snp table

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile);
GetOptions(
	   'i|infile=s'  => \$infile,
	   );


#####MAIN#####

# readin the file and store data
my $snps = {};

open (I, "$infile");
while (my $line = <I>){
    chomp $line;

    my ($acc, $seq) = 
	split (/\t/, $line);
    
    # blow apart the sequence
    my @seq = split (//, $seq);

    # plug in the snp chars
    my $pos = 0;
    foreach my $snp (@seq){
	$pos++;
	$snps->{$pos}->{$acc} = $snp;
    }
    print STDERR "$acc\t$pos\n";
}
close (I);

foreach my $pos (sort {$a <=> $b} keys %$snps){
    my $uniqbases = {};
    foreach my $acc (sort keys %{$snps->{$pos}}){
	$uniqbases->{$snps->{$pos}->{$acc}}++;
    }
    print "$pos\t";
    my @string;
    my $basecount = 0;
    foreach my $base (sort keys %$uniqbases){
	$basecount++;
	push (@string, $base, $uniqbases->{$base});
    }
    my $string = join "\t", @string;
#    print "$string\n";
    print "$basecount\n";
}

