#!/usr/bin/perl

# Purpose: Accept a file with all pw wga snp comparisons and calculate stats

# Options:
#          -i input file
#          -b is the number of partitions

# add the following libraries
#use lib "/home/apurva/perl";

# use the following modules 
use Getopt::Std;
use Statistics::Descriptive;

# get inputs                                                                    
my %opts = ();
getopts ('i:b:', \%opts);
my $infile = $opts{'i'};
my $bins   = $opts{'b'};

my @bins;
for (my $i = 50; $i <= 5000; $i += 50){
    push (@bins, $i);
}
#push (@bins, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40);

# store all the data
$data = {};
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    my ($query, $ref, $snpcount) = split (/\t/, $line);
    $data->{$query}->{$ref} = $snpcount;
}
close (F);

# average each pair
my $seen = {};
my @lens;
foreach my $q (sort keys %$data){
    foreach my $r (sort keys %{$data->{$q}}){
	next if ($q eq $r);
	my $qr = $q . $r;
	my $rq = $r . $q;
	next if ((exists ($seen->{$qr})) or (exists ($seen->{$rq})));
	if ((exists ($data->{$q}->{$r})) and (exists ($data->{$r}->{$q}))){
	    my $avg = ($data->{$q}->{$r} + $data->{$r}->{$q}) / 2;
	    my $qr = $q . $r;
	    my $rq = $r . $q;
	    $seen->{$qr} = 1;
	    $seen->{$rq} = 1;
	    push (@lens, $avg);
	}
	else {
	    print STDERR "Missing\t$q:$r\n";
	}
    }
}
# do descriptive stats
my $statobj = Statistics::Descriptive::Full->new();
$statobj->add_data(@lens);
    
my $count = $statobj->count();
my $mean  = $statobj->mean();
my $median = $statobj->median();
my $stddev = $statobj->standard_deviation();
#my %freqs = $statobj->frequency_distribution($bins);
my %freqs = $statobj->frequency_distribution(\@bins);

print "count\t$count\n";
print "mean\t$mean\n";
print "median\t$median\n";
print "stddev\t$stddev\n";

my $sum = 0;
for (sort {$a <=> $b} keys %freqs) {
    print "$_\t$freqs{$_}\n";
    $sum += $freqs{$_};
}

print "$sum\n";



    
