#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:', \%opts);
my $file = $opts{'f'};

my $counter = 0;
open (F, "$file");
while (my $line = <F>){
    $counter++;
    chomp $line;
    print "$counter\t$line\n";
#    (print "$line\n") unless ($line =~m/CHARSET/);
}
close (F);
