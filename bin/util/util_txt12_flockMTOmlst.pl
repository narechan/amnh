#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($infile, $help);
GetOptions(
    'h|help'          => \$help,
    'i|in=s'       => \$infile,
    ) or pod2usage;

pod2usage if $help;

for my $option ($infile){
    (warn ("Missing a required option\n") and pod2usage)                            
        unless ($option);                                                           
}                                                                                   

####MAIN####                                                                         

# create a symbol library (each letter of the alphabet)
my $library = {};
my $c = 0;
for my $letter ("A".."Z"){
    $c++;
    $library->{$letter} = $c
}

# parse the nexus file and store the data                                            
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    my ($acc, $data) = split (/\t/, $line);
    my @data = split (//, $data);
    print "$acc\t";
    my @string;
    foreach my $d (@data){
	push (@string, $library->{$d});
    }
    my $string = join "\t", @string;
    print "$string\n";
}
close (F);

