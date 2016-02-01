#!/usr/bin/perl -w

# -i is the input table
# -l is the length of the window you want
# -g is the length of your genome
# -m is the motion (number of bases you want to slide your window
# -r is the name of the refernce to be reflected in circos

#####SETUP#####

use strict;
use Getopt::Long;

my ($infile, $length, $genlen, $motion, $reference);
GetOptions(
    'i|infile=s'   => \$infile,
    'l|length=s'   => \$length,
    'g|genlen=s'   => \$genlen,
    'm|motion=s'   => \$motion,
    'r|reference=s'=> \$reference,
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

# cycle through the infile and tabulate
my $freqs = {};
my $counter = 0;
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    $counter++;
    print STDERR "$counter\n";

    my ($acc, $score) = 
	split (/ /, $line);
    
    # brute forcce binning
    foreach my $region (sort {$a <=> $b} keys %$parts){
	my ($start, $end) = split (/\-/, $parts->{$region});
	if (($start <= $counter) and ($counter <= $end)){
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
