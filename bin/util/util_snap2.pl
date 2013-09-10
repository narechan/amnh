#!/usr/bin/perl                                                                                                  

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:', \%opts);
my $file  = $opts{'f'};

open (F, "$file");
while (my $line = <F>){
    chomp $line;
    next if ($line =~m/^\#/);
    
    my @line = split (/\s+/, $line);
    shift @line;
    
    if ($line[5] == 0){
	print "$line[0]\tundef\n";
    }
    else{
	my $dnds = $line[6] / $line[5];
	print "$line[0]\t$dnds\n";
    }
}
