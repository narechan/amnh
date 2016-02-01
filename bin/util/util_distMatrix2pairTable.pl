#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('i:c:', \%opts);
my $infile    = $opts{'i'};
my $chars     = $opts{'c'};

my @chars;
open (C, "$chars");
while (my $line = <C>){
    chomp $line;
    push (@chars, $line);
}
close (C);

my $pairs = {};
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    my @data = split (/\t/, $line);
    my $char = shift @data;
    
    my $counter = -1;
    foreach my $val (@data){
	$counter++;
	my $mate = $chars[$counter];
	
	if (exists ($pairs->{$char}->{$mate})){
	    next;
	}
	elsif (exists ($pairs->{$mate}->{$char})){
	    next;
	}
	else {
	    $pairs->{$char}->{$mate} = $val;
	}
    }
}
close (F);

foreach my $char1 (sort keys %$pairs){
    foreach my $char2 (sort keys %{$pairs->{$char1}}){
	my $pair = $char1 . "v" . $char2;
	print "$pair\t$pairs->{$char1}->{$char2}\n";
    }
}
