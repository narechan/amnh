#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;

my ($help, $bamfile, $refseq, $outdir);
GetOptions(
    'h|help'          => \$help,
    'b|bamfile=s'     => \$bamfile,
    'r|ref=s'         => \$refseq,
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;

pod2usage if $help;

for my $option ($bamfile, $refseq, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# parse the reference
my $lengths = {};
my $totlen  = 0;
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$refseq");
while (my $sequence_obj = $seqin->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $len      = $sequence_obj->length();
    $totlen += $len;
    $lengths->{$id} = $len;
}

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
`samtools depth $bamfile > $outdir/depth.$bam_name.out`;

# parse the depth file
my $sum = {};
my $touched = {};
my $sumtot = 0;
my $touchtot = 0;
open (D, "$outdir/depth.$bam_name.out");
while (my $line = <D>){
    chomp $line;
    my ($acc, $pos, $cov) = split (/\t/, $line);
    $sum->{$acc} += $cov;
    $touched->{$acc}++;
    $sumtot += $cov;
    $touchtot++;
}
close (D);

open (OUT, ">$outdir/$bam_name.coverage");
foreach my $contig (keys %$lengths){
    my $coverage;
    my $touchy;
    if ($sum->{$contig}){
	$coverage = $sum->{$contig} / $lengths->{$contig};
    }
    else {
	$coverage = 0;
    }
    if ($touched->{$contig}){
	$touchy  = $touched->{$contig} / $lengths->{$contig};
    }
    else{
	$touchy = 0;
    }
    print OUT "$contig\t$lengths->{$contig}\t$coverage\t$touchy\n";
}

my $c = $sumtot / $totlen;
my $t = $touchtot / $totlen;
print OUT "$totlen\t$c\t$t\n";
close (OUT);
