#!/usr/bin/perl

# Purpose: Accept a directory of jk summaries and generate stats

# Options:
#          -i input dir
#          -l is the list of things you want to analyze


# use the following modules 
use Getopt::Std;
use Statistics::Descriptive;

# get inputs                                                                    
my %opts = ();
getopts ('i:l:', \%opts);
my $indir = $opts{'i'};
my $list  = $opts{'l'};

# get the list
my %list;
open (F, "$list");
while (my $line = <F>){
    chomp $line;
    $list{$line} = 1;
}
close (F);

# get the files
opendir (E, "$indir");
my @sums = sort (readdir(E));
shift @sums;
shift @sums;
closedir (E);

# parse the data
my %perexcl;
foreach my $sum (@sums){
    next if ($sum =~m/stderr/);
    next unless (exists ($list{$sum}));
    warn "$sum\n";
    open (F, "$indir/$sum/summary");
    my ($rep, $crap) = split (/\./, $sum);
    while (my $line = <F>){
	chomp $line;
	my ($excl, $scstr, $avesc) = split (/\t/, $line);
	push @{$perexcl{$excl}}, $avesc;
	print "$rep\t$line\n";
    }
    close (F);
}

foreach my $excl (sort {$a <=> $b} keys %perexcl){
    my $statobj = Statistics::Descriptive::Full->new();
    $statobj->add_data(@{$perexcl{$excl}});
    
    my $mean  = $statobj->mean();
    my $stddev = $statobj->standard_deviation();
    
    print STDERR "$mean\t$stddev\n";
}
