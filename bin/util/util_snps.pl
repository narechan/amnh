#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile);
GetOptions(
	   'i|infile=s'  => \$infile,
	   );

#####MAIN#####

my $vrsa = 0;
my $vssa = 0;
my $both = 0;
my $n = 0;
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    next if ($line =~m/^Contig/);
    my ($refctg, $refpos, $vrsab, $vssab, $refb, $state, $stuff) = 
	split (/\t/, $line, 7);

    if (($vrsab eq $vssab) and ($vrsab ne $refb) and ($vssab ne $refb)){
	$n++;
    }
    if ($vrsab ne $refb){
	$vrsa++;
    }
    if ($vssab ne $refb){
        $vssa++;
    }
    if ($vssab ne $vrsab){
	$both++
    }
}
print "$vrsa\t$vssa\t$both\t$n\n";
close (I);

   
