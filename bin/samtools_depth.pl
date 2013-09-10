#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $bamfile, $size);
GetOptions(
    'h|help'          => \$help,
    'b|bamfile=s'     => \$bamfile,
    's|size=i'         => \$size,	   
    ) or pod2usage;

pod2usage if $help;

for my $option ($bamfile, $size){
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
`samtools depth $bamfile > /tmp/depth.$bam_name.out`;

# parse the depth file
my $sum = 0;
open (D, "/tmp/depth.$bam_name.out");
while (my $line = <D>){
    chomp $line;
    my ($acc, $pos, $cov) = split (/\t/, $line);
    $sum += $cov;
}
close (D);

my $coverage = $sum / $size;
print "$bam_name\t$coverage\n";
`rm /tmp/depth.$bam_name.out`;
