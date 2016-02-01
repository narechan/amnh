#!/usr/bin/perl -w

# Generate an asymptotic snp collector curve as taxa are added and compared to a single reference

# -i is an input directory of query genomes
# -r is the reference genome
# -o is the output dir

#####SETUP#####

use strict;
use Getopt::Long;

my ($indir, $ref, $outdir);
GetOptions(
	   'i|indir=s'   => \$indir,
	   'r|ref=s'   => \$ref,
	   'o|outdir=s'   => \$outdir,
	   );

`mkdir -p $outdir`;

#####MAIN#####

# load all the query genomes
opendir (D, "$indir");
my @querygenomes = readdir (D);
shift @querygenomes;
shift @querygenomes;
closedir (D);

# cycle through the query genomes randomly tracking snp positions
my $reftracker = {};
while (@querygenomes){

    # randomly select something
    my $randindex = rand (@querygenomes);
    my $randgenome = $querygenomes[$randindex];
    
    # align and find the snps
    `nucmer --prefix=$outdir/$randgenome $ref $indir/$randgenome`;
    `show-snps -CTHr $outdir/$randgenome.delta > $outdir/$randgenome.snps`;

    # parse the snps out
    my $snps = 0;
    open (S, "$outdir/$randgenome.snps");
    while (my $line = <S>){
	chomp $line;
	my @line = split (/\t/, $line);
	next if ($line[1] eq ".");
	next if ($line[2] eq ".");
	$snps++;
	$reftracker->{$line[8]}->{$line[0]}->{$line[2]} = 1;
    }
    close (S);

    # output the snp increments
    my $refposchanges = 0;
    my $totalalleles = 0;
    foreach my $contig (sort keys %$reftracker){
	foreach my $pos (sort keys %{$reftracker->{$contig}}){
	    $refposchanges++;
	    foreach my $allele (sort keys %{$reftracker->{$contig}->{$pos}}){
		$totalalleles++;
	    }
	}
    }
    print "$randgenome\t$snps\t$refposchanges\t$totalalleles\n";

    # remove the randomly selected genome
    splice (@querygenomes, $randindex, 1);
}
	

