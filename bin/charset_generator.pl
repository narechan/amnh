#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:n:', \%opts);
my $interval  = $opts{'i'};
my $nchars    = $opts{'n'};

my $start = 1;
my $end   = $interval;
print "BEGIN SETS;\n";
my $counter = 0;
until ($start > $nchars){
    $counter++;
    my $end = $start + $interval - 1;
    
    if ($end > $nchars){
	print "CHARSET C$counter = $start - $nchars;\n";
    }
    else {
	print "CHARSET C$counter = $start - $end;\n";
    }

    $start = $end + 1;
}
print "END;\n";
