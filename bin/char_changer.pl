#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('f:', \%opts);
my $file = $opts{'f'};

my $counter = 0;
open (F, "$file");
while (my $line = <F>){
    chomp $line;
    ($line =~s/\[/\n\[/) if ($line =~m/^\[/);
    print "$line\n";
}
close (F);
