#!/usr/bin/perl

# Purpose: Accept a directory of ls summaries and generate stats

# Options:
#          -i input dir
#          -o outdir


# use the following modules 
use Getopt::Std;
use Statistics::Descriptive;

# get inputs                                                                    
my %opts = ();
getopts ('i:o:', \%opts);
my $indir = $opts{'i'};
my $outdir = $opts{'o'};

`mkdir -p $outdir`;

# get the files
opendir (E, "$indir");
my @sums = sort (readdir(E));
shift @sums;
shift @sums;
closedir (E);

# parse the data
my $data = {};
open (A, ">$outdir/all");
foreach my $sum (@sums){
    my ($rep, $parts, $stuff) = split (/\./, $sum, 3);
    
    open (F, "$indir/$sum/summary");
    while (my $line = <F>){
	chomp $line;
	my ($node, $neg, $pos, $diff, $stat) = split (/\t/, $line);
	
	push @{$data->{$node}->{$parts}}, $stat;
	print A "$rep\t$parts\t$line\n";
    }
    close (F);
}
close (A);

# do descriptive stats
open (S, ">$outdir/sum");
my $counter = 0;
foreach my $node (sort {$a <=> $b} keys %$data){
    foreach my $part (sort {$a <=> $b} keys %{$data->{$node}}){

	my $statobj = Statistics::Descriptive::Full->new();
	$statobj->add_data(@{$data->{$node}->{$part}});
    
	my $count = $statobj->count();
	my $mean  = $statobj->mean();
	
	print S "$node\t$part\t$count\t$mean\n";
    }
}	
close (S);
