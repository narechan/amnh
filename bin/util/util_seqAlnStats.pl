#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:d:s:', \%opts);
my $infile = $opts{'i'};
my $depthdir  = $opts{'d'};
my $snpoutdir = $opts{'s'};

my $denom = 4411532;

open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    my ($strain, $reads) = split (/\t/, $line);
    
    open (D, "$depthdir/$strain-h37rv.depth");
    my $counter = 0;
    my $hittot = 0;
    while (my $line = <D>){
	chomp $line;
	$counter++;
	
	my ($ref, $pos, $hits) = split (/\t/, $line);
	$hittot += $hits;
    }
    close (D);
    
    my $readdepth = $hittot / $counter;
    my $coverage  = $counter / $denom;
    
    open (S, "$snpoutdir/$strain-h37rv.snpout");
    my $snpcounter = 0;
    my $filtsnps = 0;
    while (my $line = <S>){
	chomp $line;
	$snpcounter++;
	
	my @line = split (/\t/, $line);
	($filtsnps++) if ($line[5] >= 20);
    }
    close (S);
    
    print "$line\t$readdepth\t$coverage\t$snpcounter\t$filtsnps\n";
}
close (I);
	  
