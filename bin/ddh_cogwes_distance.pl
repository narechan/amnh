#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('f:t:', \%opts);
my $file  = $opts{'f'};
my $taxa  = $opts{'t'};

#my $refdata = {};
my $qrydata = {};
open (F, "$file");
while (my $line = <F>){
    chomp $line;
    my @data = split (/\t/, $line);
    
#    next unless ($data[2] eq "solexa");

    $qrydata->{$data[0]}->{$data[1]} = $data[9];
#    $refdata->{$data[1]}->{$data[0]} = $data[9];
}
close (F);

my @srtqrydata = sort (keys %$qrydata);

foreach my $qry (@srtqrydata){
    next unless ($qry =~m/$taxa/);
    print "$qry\t";
    
    foreach my $ref (@srtqrydata){
	next unless ($ref =~m/$taxa/);
	if ($ref eq $qry){
	    print "0\t";
	}
	else{
	    my $qrycog = $qrydata->{$qry}->{$ref};
	    my $refcog = $qrydata->{$ref}->{$qry};
	    
	    my $avgcog = ($qrycog + $refcog) / 2;
	    my $cogdist = 1 - $avgcog;
	    
	    print "$cogdist\t";
	}
    }
    print "\n";
}
