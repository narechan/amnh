#!/usr/bin/perl -w

# This program accepts a list of snps from snp_annotator.pl
# and creates a density flatfile that gives the freq of snps
# in a user input interval.

# -i is the input snp table
# -l is the length of the window you want
# -t is the list of taxa you want to include
#    (others in the snp table will be skipped)
# -q is the snp quality above which snps are accepted
# -g is the length of your genome
# -m is the motion (number of bases you want to slide your window

#####SETUP#####

use strict;
use Getopt::Long;

my ($infile, $taxalist, $quality, $length, $genlen, $motion);
GetOptions(
	   'i|infile=s'   => \$infile,
	   't|taxalist=s' => \$taxalist,
	   'q|quality=s'  => \$quality,
	   'l|length=s'   => \$length,
	   'g|genlen=s'   => \$genlen,
	   'm|motion=s'   => \$motion,
	   );

#####MAIN#####

# readin the quey list
my $queries = {};
open (L, "$taxalist");
while (my $line = <L>){
    chomp $line;
    $queries->{$line} = 1;
}
close (L);

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
open (I, "$infile");
#my $scount = 0;
my $reference;
while (my $line = <I>){
    chomp $line;
#    $scount++;
#    print STDERR "$scount\n";

    my ($ref, $query, $refpos, $refbase, $snp, $sq, $stuff) = 
	split (/\t/, $line, 7);
    $reference = $ref;

    # filters
    next if ($refbase eq "N"); #ambiguous ref
    next unless (exists ($queries->{$query})); # unincluded query
    if ($sq < $quality){ #sq
	next;
    }
    if ($snp =~m/\,/){ #hets
	next;
    }

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
