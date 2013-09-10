#!/usr/bin/perl

# Purpose: Accept a directory of ts summaries and generate stats

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

# get the files
opendir (E, "$indir");
my @sums = sort (readdir(E));
shift @sums;
shift @sums;
closedir (E);

# parse the data
open (A, ">$outdir/all");
my @bs;
my @pbs;
my @phbs;
foreach my $sum (@sums){
    next if ($sum =~m/stderr/);

    open (F, "$indir/$sum/summary");
    while (my $line = <F>){
	next if ($line =~m/^Part/);
	chomp $line;
	my @data = split (/\t/, $line);
	unless ($data[3] eq "skipped"){
	    print A "$data[9]\t";
	    print A "$data[13]\t";
	    print A "$data[15]\n";
	    push (@bs, $data[9]);
	    push (@pbs, $data[13]);
	    push (@phbs, $data[15]);
	}
    }
    close (F);
}
close (A);


# do descriptive stats
my $start = -300;
my $end = 300;
my @bins;
until ($start > $end){
    push (@bins, $start);
    $start += 10;
}

open (S, ">$outdir/sum");
my $counter = 0;
for my $stat (\@bs, \@pbs, \@phbs){
    $counter++;
    my $statobj = Statistics::Descriptive::Full->new();
    $statobj->add_data(@$stat);
    
    my %freqs = $statobj->frequency_distribution(\@bins);
    my $count = $statobj->count();
    my $mean  = $statobj->mean();
    my $max = $statobj->max();
    my $min = $statobj->min();
	
    foreach my $bin (sort {$a <=> $b} keys %freqs) {
	my $ratio = $freqs{$bin} / $count;
	print S "$bin\t$freqs{$bin}\t$ratio\n";
    }
    
    print S "$mean\t";
    print S "$max\t";
    print S "$min\t";
    print S "$count\n";
}	
close (S);
