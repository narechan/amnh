#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $bamfile, $size, $window, $chrom);
GetOptions(
    'h|help'          => \$help,
    'b|bamfile=s'     => \$bamfile,
    's|size=i'         => \$size,
    'w|window=s'    => \$window,
    'c|chrom=s'     => \$chrom,
    ) or pod2usage;

pod2usage if $help;

for my $option ($bamfile, $size, $window, $chrom){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# get the bamfile name
my $bam_name;
if ($bamfile =~/\//g){
    $bamfile =~m/.*\/(.*)$/;
    $bam_name = $1;
}

else {
    $bam_name = $bamfile;
}

# do the samtools depth calculation
`samtools depth $bamfile > /scratch/apurva/tmp/depth.$bam_name.out`;

# parse the depth file
my $freq = {};
open (D, "/scratch/apurva/tmp/depth.$bam_name.out");
while (my $line = <D>){
    chomp $line;
    my ($acc, $pos, $cov) = split (/\t/, $line);
    next unless ($acc eq $chrom);
    $freq->{$pos} = $cov;
}
close (D);


my $subcounter = 0;
my $sum = 0;
for (my $i = 1; $i <= $size; $i++){
    $subcounter++;
    
    if ($subcounter == $window){
	my $coverage = $sum / $subcounter;
	my $start = $i - $window + 1;
	print "$chrom\t$start\t$i\t$coverage\n";
	$subcounter = 0;
	$sum = 0;
    }
    elsif ($i == $size){
	my $coverage = $sum / $subcounter;
	my $start = $i - $window + 1;
	print "$chrom\t$start\t$i\t$coverage\n";
        $subcounter = 0;
        $sum = 0;
    }
    else{
	if (exists ($freq->{$i})){
	    $sum += $freq->{$i};
	}
	else{
	    next;
	}
    }
}
close (D);

`rm /scratch/apurva/tmp/depth.$bam_name.out`;
