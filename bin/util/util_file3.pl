#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('f:', \%opts);
my $file  = $opts{'f'};

open (F, "$file");
while (my $line = <F>){
    chomp $line;
    if ($line =~m/\=/){
	print "CHARSET $line\n";
    }
    else {
	print "$line\n";
    }
}
close (F);


