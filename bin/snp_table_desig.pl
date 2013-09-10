#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:', \%opts);
my $in    = $opts{'i'};

open (I, "$in");
while (my $line = <I>){
    chomp $line;
    if ($line =~m/Contig\tSNP/){
	print "$line\n";
	next;
    }

    my @data = split (/\t/, $line);
    my @data2 = @data;
    pop @data2;
    shift @data2;
    shift @data2;
    
    my $rate;
#    if ($data[2] eq $data[3] eq $data[4] eq $data[5]){
#    if ($data[2] eq $data[3] eq $data[4]){
    if (keys %{{ map {$_, 1} @data2 }} == 1) {
	$rate = "N";
    }
    else {
	if (($data[2] eq "?") or ($data[3] eq "?") or ($data[4] eq "?") or ($data[5] eq "?")){
	    $rate = "+";
	}
	else {
	    $rate = "*";
	}
    }
    
    print "$line\t$rate\n";
    
}
close (I);
