#!/usr/bin/perl -w

=head1 NAME

bfastSamtools_runner.pl

=head1 SYNOPSIS

  bfastSamtools_runner.pl -- this program executes bwa and samtools for variant calling

Options:

 --help        Show brief help and exit
 --index       Contains the bfast index
 --seq         Contains your reads
 --name        Is the prefix for your output files
 --cs          Indicates that you need to aln in colorspace
 --outdir      Is your output dir

bfast and samtools must be in the path

changes to the parameters to either bwa or samtools must be entered 
into the script below (HARD CODED)

=head1 DESCRIPTION

Run bfast and samtools to find variants.

Usage examp:

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####
#####     

#####SETUP#####
    
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $index, $seq, $outdir, $name, $cs);
GetOptions(
    'h|help'          => \$help,
    'i|index=s'       => \$index,
    's|seq=s'           => \$seq,
    'o|outdir=s'      => \$outdir,
    'n|name=s'        => \$name,
    'c|cs'          => \$cs,
    ) or pod2usage;

pod2usage if $help;

for my $option ($index, $seq, $outdir, $name){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/results`;
`mkdir -p $outdir/stderr`;

#####MAIN#####

# get filenames
my $indexname = find_name ($index);
my $seqname   = find_name ($seq);

# align reads to the reference
if ($cs){
    `bfast match -l -f $index -A 1 -n 8 -r $seq 1>$outdir/results/$name.bmf 2>$outdir/stderr/$name.match`;
    `bfast localalign -n 8 -f $index -A 1 -m $outdir/results/$name.bmf 1>$outdir/results/$name.baf 2>$outdir/stderr/$name.localalign`;
    `bfast postprocess -f $index -i $outdir/results/$name.baf -A 1 1>$outdir/results/$name.sam 2>$outdir/stderr/$name.pp`;
}
else{
    `bfast match -l -f $index -n 8 -r $seq 1>$outdir/results/$name.bmf 2>$outdir/stderr/$name.match`;
    `bfast localalign -n 8 -f $index -m $outdir/results/$name.bmf 1>$outdir/results/$name.baf 2>$outdir/stderr/$name.localalign`;
    `bfast postprocess -f $index -i $outdir/results/$name.baf 1>$outdir/results/$name.sam 2>$outdir/stderr/$name.pp`;
}
    
# convert to bam and sort
`samtools view -bS $outdir/results/$name.sam | samtools sort - $outdir/results/$name.srt 2>$outdir/stderr/$name.bam`;

# Index bam file
`samtools index $outdir/results/$name.srt.bam 2>$outdir/stderr/$name.index`;

# call raw snps -- note that -f is referring to a fasta file of the reference
# so that should be in there too.
#`samtools pileup -cv -f $index $outdir/results/$name.srt.bam 1>$outdir/results/$name.pileup 2>$outdir/stderr/$name.pileup`;
`samtools mpileup -uf $index $outdir/results/$name.srt.bam | bcftools view -bvcg - 1>$outdir/results/$name.bcf 2>$outdir/stderr/$name.bcf`;
`bcftools view $outdir/results/$name.bcf > $outdir/results/$name.vcf`;

# parse the vcf
open (PU, "$outdir/results/$name.vcf");
open (MI, ">$outdir/results/$name.snpsonly");
open (ID, ">$outdir/results/$name.indels");
while (my $line = <PU>){
    next if ($line =~m/^\#/);
    chomp $line;
    my ($refseq, $refpos, $rsid, $refbase, $snp, $sq, $filter, $info, $format, $stuff)
	= split (/\t/, $line);
    if ($info =~m/INDEL/){
	print ID "$line\n";
    }
    else{
	print MI "$line\n";
    }
}
close (ID);
close (MI);
close (PU);

 

#####SUBS#####

sub find_name{
    my $name = shift;
    
    my $fasta_name;
    if ($name =~/\//g){
	$name =~m/.*\/(.*)$/;
	$fasta_name = $1;
    }
    
    else {
	$fasta_name = $name;
    }
    return ($fasta_name);
}
