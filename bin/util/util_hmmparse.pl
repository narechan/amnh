#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:', \%opts);
my $infile  = $opts{'i'};

my $a = {};
my $b = {};
open (A, "$infile");
while (my $line = <A>){
    chomp $line;
    next if ($line =~m/\#/);
    my @line = split (/\s+/, $line);
    
    if (exists ($a->{$line[2]})){
	if ($line[5] > $b->{$line[2]}){
	    $a->{$line[2]} = $line[0];
	    $b->{$line[2]} = $line[5];
	}
	else {
	    next;
	}
    }
    else {
	$a->{$line[2]} = $line[0];
	$b->{$line[2]} = $line[5];
    }
}
close (A);

foreach my $pt (keys %$a){
    print "$pt\t$a->{$pt}\n";
#    my $ptind =  $pt;
#    $ptind =~s/SA//g;
    
#    if ((($ptind >= 2339) and ($ptind <=2502
}

