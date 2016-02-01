#!/usr/bin/perl

use Getopt::Std;
use PDL::IO::HDF5;

my %opts = ();
getopts ('i:', \%opts);
my $infile = $opts{'i'};

my $hdf = new PDL::IO::HDF5("$infile");
my @groups = $hdf->groups;
my $group = $hdf->group("AlnInfo");
my @attrs = $group->attrs;
my @datasets = $group->datasets;
my $dataset = $group->dataset("AlnIndex");
my @dims = $dataset->dims;
my $pdl = $dataset->get(); 
my $n = $pdl->nelem();
my $info = $pdl->info;
my $num = $pdl->at(5,3);
my $slice_1b = $pdl->slice("0:2,-1");
print "$slice_1b\n";
print yellow;
