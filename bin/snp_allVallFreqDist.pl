#!/usr/bin/perl -w

# do all versus all pw wga for distribution of snps frequencies for a set of genomes

# -i is an input directory of query genomes
# -r is the reference
# -o is the output dir

#####SETUP#####

use strict;
use Getopt::Long;

my ($indir, $outdir, $ref);
GetOptions(
	   'i|indir=s'   => \$indir,
	   'o|outdir=s'   => \$outdir,
           'r|ref=s'     => \$ref,
	   );

`mkdir -p $outdir`;

#####MAIN#####

# load all the query genomes
opendir (D, "$indir");
my @querygenomes = readdir (D);
shift @querygenomes;
shift @querygenomes;
closedir (D);

my $refname;
if ($ref =~/\//g){
    $ref =~m/.*\/(.*)$/;
    $refname = $1;
}

else {
    $refname = $ref;
}

# cycle through the query genomes and do pw wga aligns
foreach my $q1 (@querygenomes){
    
    # align and find the snps
    `nucmer --prefix=$outdir/$q1-$refname $ref $indir/$q1`;
    `show-snps -CTHr $outdir/$q1-$refname.delta > $outdir/$q1-$refname.snps`;

    # parse the snps out
    my $snps = 0;
    open (S, "$outdir/$q1-$refname.snps");
    while (my $line = <S>){
	chomp $line;
	my @line = split (/\t/, $line);
	next if ($line[1] eq ".");
	next if ($line[2] eq ".");
	$snps++;
    }
    close (S);

    print "$q1\t$refname\t$snps\n";
    `rm $outdir/*`;
}
	

