#!/usr/bin/perl -w

# This program accepts a list of snps from snpeff
# and creates a density flatfile that gives the freq of snps
# in a user input interval.

# -i is the input snp table or vcf file (a single snpeff or vcf, not concatenated across samples)
#    if you want to compile stats across queries you can concatenate them but their origin will be lost
# -l is the length of the window you want
# -g is the length of your genome
# -m is the motion (number of bases you want to slide your window

# NOTE that snp calls should already be at the filtered stage here (hets, cutoff, omitted regions, etc...)
# you can do that in snp_matrixBuilder.pl

# either -i or -v must be specified

#####SETUP#####

use strict;
use Getopt::Long;

my ($infile, $length, $genlen, $motion, $vcffile);
GetOptions(
	   'i|infile=s'   => \$infile,
	   'l|length=s'   => \$length,
	   'g|genlen=s'   => \$genlen,
	   'm|motion=s'   => \$motion,
           'v|vcffile=s'   => \$vcffile,
	   );

#####MAIN#####

# generate windows
my $parts = {};
my $counter = 0;
my $i = 1;
while (1){
    $counter++;
    my $start = $i;
    my $end   = $i + ($length - 1);

    my $coords = $start . "-" . $end;

    if ($end >= $genlen){
	$end = $genlen;
	$coords = $start . "-" . $end;
	$parts->{$counter} = $coords;
	last;
    }

    $parts->{$counter} = $coords;
    $i += $motion;
}

# cycle through the snps file and tabulate
my $freqs = {};
my $reference;

open (I, "$infile");
#my $scount = 0;
while (my $line = <I>){
    chomp $line;
#    $scount++;
#    print STDERR "$scount\n";
    
    my ($ref, $refpos, $stuff) = 
	split (/\t/, $line, 3);
    $reference = $ref;
    
    # brute forcce binning
    foreach my $region (sort {$a <=> $b} keys %$parts){
	my ($start, $end) = split (/\-/, $parts->{$region});
	if (($start <= $refpos) and ($refpos <= $end)){
	    $freqs->{$region}++;
	    last;
	}
	else {
	    next;
	}
    }
}
close (I);
    
# fill in all undefined values in freq datastruc
foreach my $region (sort {$a <=> $b} keys %$parts){
    if (exists ($freqs->{$region})){
	next;
    }
    else {
	$freqs->{$region} = 0;
    }
}

# print the freq table
my $total = 0;
foreach my $bin (sort {$a <=> $b} keys %$freqs){
    my $region = $parts->{$bin};
    my ($start, $end) = split (/\-/, $region);
    
    print "$reference\t$start\t$end\t$freqs->{$bin}\n";
    $total += $freqs->{$bin};
}
print STDERR "$total\n";
