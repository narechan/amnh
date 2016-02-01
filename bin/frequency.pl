#!/usr/bin/perl

# Purpose: Accept a single column data file and generate sequence stats

# Options:
#          -i input file
#          -b is the number of partitions

# add the following libraries
use lib "/home/apurva/perl";

# use the following modules 
use Getopt::Std;
use Statistics::Descriptive;
use Bio::Seq;
use Bio::SeqIO;

# get inputs                                                                    
my %opts = ();
getopts ('i:b:', \%opts);
my $infile = $opts{'i'};
my $bins   = $opts{'b'};

my @lens;
my @bins;
push (@bins, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40);
#my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$infile");
#while (my $sequence_obj = $seqin->next_seq()){
#    my $id       = $sequence_obj->display_id();
#    my $seq      = $sequence_obj->seq();
#    my $len      = length ($seq);
#    push (@lens, $len);
#}
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    push (@lens, $line);
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



    
