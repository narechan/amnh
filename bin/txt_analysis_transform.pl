#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('f:', \%opts);
my $file  = $opts{'f'};

my $text;
open (F, "$file");
while (my $line = <F>){
    chomp $line;
    $line =~tr/[a-z]/[A-Z]/;

    my @line = split (//, $line);
    foreach my $char (@line){
	next unless ($char =~m/[a-zA-Z]/);
	if ($char =~m/B/){
	    print "C";
	}
	elsif ($char =~m/J/){
            print "K";
	}
	elsif ($char =~m/O/){
	    print "P";
	}
	elsif ($char =~m/U/){
	    print "V";
	}
	elsif ($char =~m/X/){
	    print "Y";
	}
	elsif ($char =~m/Z/){
	    print "A";
	}
	else{
	    print $char;
	}
    }
}
print "\n";
close (F);
