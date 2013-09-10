#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('f:', \%opts);
my $file  = $opts{'f'};

open (F, "$file");
while (my $line = <F>){
    chomp $line;
    my ($parts, $walltime) = split (/\t/, $line);
    my ($minutes, $seconds) = split (/m/, $walltime);
    $seconds =~s/s//g;
    
    my $minfrac = $seconds / 60;
    my $totmins = $minutes + $minfrac;
    
    my $cpuhours = $totmins / 60;
    
    print "$parts\t$walltime\t$cpuhours\n";
}
close (F);

    
