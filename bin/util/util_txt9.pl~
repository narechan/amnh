#!/usr/bin/perl

use Getopt::Std;

my %opts = ();
getopts ('l:i:', \%opts);
my $list  = $opts{'l'};
my $indir    = $opts{'i'};

opendir (D, "$indir");
my @exptfiles = sort (readdir(D));
shift @exptfiles;
shift @exptfiles;
closedir (E);

my $listdata = {};
open (F, "$list");
while (my $line = <F>){
    chomp $line;
    my @data = split (/\t/, $line);
    my @three = split (/\//, $data[2]);
    my $val = pop @three;
    my $val2 = pop @three;
    $listdata->{$val2} = $data[1];
}
close (F);

foreach my $file (@exptfiles){
    my ($prefix, $junk) = split (/\./, $file);
    `mv $indir/$file $indir/$listdata->{$prefix}.fastq`;
}
