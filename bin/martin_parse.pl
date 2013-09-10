#!/usr/bin/perl

use Getopt::Std;
use Array::Transpose;
use Data::Dumper;

my %opts = ();
getopts ('f:h:', \%opts);
my $file  = $opts{'f'};
my $headers = $opts{'h'};

$/ = "\n\n";

my $total = {};
open (F, "$file");
while (my $chunk = <F>){
    chomp $chunk;
    next if ($chunk =~m/PROSOFT_INSTRUMENT_CONFIG/);

    # build chunk matrix
    my @master;
    my @data = split (/\n/, $chunk);
    my $counter = -1;
    foreach my $row (@data){
	$counter++;
	my @row = split ("\t", $row);
	push @master, [@row];
    }
    
    #transpose the matrix
    my @mastert = transpose (\@master);

    # populate the file datastruc
    foreach $array (@mastert){
	my $a = shift @$array;
	my $b = shift @$array;

	$total->{$a}->{$b} = $array;
    }
}
close (F);

$/ = "\n";
open (H, "$headers");
while (my $line = <H>){
    chomp $line;
    
    foreach my $b (sort (keys %{$total->{$line}})){
	print "$line\t$b\t";
	my $elestring = join ("\t", @{$total->{$line}->{$b}});
	print "$elestring\n";
    }
}
close (H);

