#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:', \%opts);
my $infile    = $opts{'i'};

my $table = {};
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    my @data = split (/\t/, $line);
    ($table->{$data[1]}++) if ($data[5] >= 20);
}
close (F);

foreach my $ent (sort keys %$table){
    print "$ent\t$table->{$ent}\n";
}
