#!/usr/bin/perl -w

# -i is the input table
# -l is the length of the window you want
# -g is the length of your genome
# -r is the name of the refernce to be reflected in circos

#####SETUP#####

use strict;
use Getopt::Long;

my ($infile, $length, $genlen, $reference);
GetOptions(
    'i|infile=s'   => \$infile,
    'l|length=s'   => \$length,
    'g|genlen=s'   => \$genlen,
    'r|reference=s'=> \$reference,
	   );

#####MAIN#####

# cycle through the infile and tabulate
my $freqs = {};
my $counter = 0;
my $mastercounter = 0;
my $sum = 0;
my $start = 1;
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    $counter++;
    $mastercounter++;
#    print STDERR "$mastercounter\n";

    my ($acc, $score) = 
	split (/\t/, $line);
    
    if ($counter == $length){
	($sum += $score) unless ($score eq "NA");
	my $coords = $start . "-" . $mastercounter;
	my $avg    = $sum / $length;
	$freqs->{$coords} = $avg;
	    
	$sum = 0;
        $counter = 0;
	$start = $mastercounter + 1;
    }
    elsif ($mastercounter == $genlen){
	($sum += $score) unless ($score eq "NA");
        my $coords = $start . "-" . $mastercounter;
        my $avg    = $sum / ($mastercounter - $start + 1);
        $freqs->{$coords} = $avg;

        $sum = 0;
        $counter = 0;
        $start = $mastercounter + 1;
    }
    else {
	($sum += $score) unless ($score eq "NA");
    }
}
close (I);

# print the freq table
my $total = 0;
foreach my $bin (sort keys %$freqs){
    my ($start, $end) = split (/\-/, $bin);
    
    print "$reference\t$start\t$end\t$freqs->{$bin}\n";
}