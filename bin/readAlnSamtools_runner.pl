#!/usr/bin/perl -w

=head1 NAME

readAlnSamtools_runner.pl

=head1 SYNOPSIS

  readAlnSamtools_runner.pl -- this program executes bwa and samtools for variant calling

Options:

 --help        Show brief help and exit
 --index       Contains the bwa index (also name of the fasta reference)
 --seq1        Contains the first set of paired ends or the single end reads (fastq)
 --seq2        Contains the second set of paired ends (fastq) (optional)
 --name        Is the prefix for your output files
 --outdir      Is your output dir
 --config      Is your config file that specifies bwa, bowtie2, or soap and optional params
 
your short read aligner and samtools must be in the path

changes to the parameters for samtools must be entered 
into the script below (HARD CODED)

NOTE: if using soap, SQ headers are not generated, so
you must provide a file that with each reference name and
its length (i.e. *.fai file in bwa index; -t option to samtools view)

=head1 DESCRIPTION

Run a short read aligner and samtools to find variants.

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

my ($help, $index, $seq1, $seq2, $outdir, $name, $snps, $configfile);
GetOptions(
    'h|help'          => \$help,
    'i|index=s'       => \$index,
    'seq1|a=s'           => \$seq1,
    'seq2|b=s'           => \$seq2,
    'o|outdir=s'      => \$outdir,
    'n|name=s'        => \$name,
    'c|configfile=s'          => \$configfile,
    ) or pod2usage;

pod2usage if $help;

for my $option ($index, $seq1, $outdir, $name, $configfile){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/results`;
`mkdir -p $outdir/stderr`;

#####MAIN#####

# get the alignner and config
my $conf = parse_config ($configfile);

# align reads to the reference
if (exists ($conf->{'BWA'})){
    if ($seq2){
	`bwa aln $conf->{'BWA'} $index $seq1 1>$outdir/results/$name-1.sai 2>$outdir/stderr/$name-1.aln`;
	`bwa aln $conf->{'BWA'} $index $seq2 1>$outdir/results/$name-2.sai 2>$outdir/stderr/$name-2.aln`;
	`bwa sampe $index $outdir/results/$name-1.sai $outdir/results/$name-2.sai $seq1 $seq2 1>$outdir/results/$name.sam 2>$outdir/stderr/$name.sam`;
    }
    else{
	`bwa aln $conf->{'BWA'} $index $seq1 1>$outdir/results/$name-1.sai 2>$outdir/stderr/$name-1.aln`;
	`bwa samse $index $outdir/results/$name-1.sai $seq1 1>$outdir/results/$name.sam 2>$outdir/stderr/$name.sam`;
    }
}

elsif (exists ($conf->{'BOWTIE2'})){
    if ($seq2){
	`bowtie2 $conf->{'BOWTIE2'} -x $index -1 $seq1 -2 $seq2 1>$outdir/results/$name.sam 2>$outdir/stderr/$name.sam`;
    }
    else {
	`bowtie2 $conf->{'BOWTIE2'} -x $index -U $seq1 1>$outdir/results/$name.sam 2>$outdir/stderr/$name.sam`;
    }
}

elsif (exists ($conf->{'SOAP'})){
    if ($seq2){
	`soap $conf->{'SOAP'} -a $seq1 -b $seq2 -D $index -o $outdir/results/$name.soap -2 $outdir/results/$name.unpaired 2>$outdir/stderr/$name.soap`;
	`soap2sam.pl -p $outdir/results/$name.soap 1>$outdir/results/$name.sam 2>$outdir/stderr/$name.sam`;
    }
    else {
	`soap $conf->{'SOAP'} -a $seq1 -D $index -o $outdir/results/$name.soap -2 $outdir/results/$name.unpaired 2>$outdir/stderr/$name.soap`;
	`soap2sam.pl -p $outdir/results/$name.soap 1>$outdir/results/$name.sam 2>$outdir/stderr/$name.sam`;
    }
    
    # bam sort specific here since soap does not make sam files with headers
    `samtools view -bSt $conf->{'INDEX'} $outdir/results/$name.sam | samtools sort - $outdir/results/$name.srt 2>$outdir/stderr/$name.bam`;
}

else {
    print STDERR "Unknown alignment method\n";
    die;
}


# convert to bam and sort
unless (exists ($conf->{'SOAP'})){
    `samtools view -bS $outdir/results/$name.sam | samtools sort - $outdir/results/$name.srt 2>$outdir/stderr/$name.bam`;
}

# index bam file
`samtools index $outdir/results/$name.srt.bam 2>$outdir/stderr/$name.index`;

# call raw snps
#`samtools pileup -cv -f $reference $outdir/results/$name.srt.bam 1>$outdir/results/$name.pileup 2>$outdir/stderr/$name.pileup`;
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


####SUB####

sub parse_config{
    my $file = shift;

    my %config;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }
    
    return (\%config);
}
