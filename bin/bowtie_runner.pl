#!/usr/bin/perl -w

=head1 NAME

bowtie_runner.pl

=head1 SYNOPSIS

  bowtie_runner.pl -- this program will batch jobs to raxml

Options:

 --help        Show brief help and exit
 --index       Contains the bowtie index
 --sequence    Contains the file of sequence
    (if use something other than fastq, specify in config options to bowtie)
 --reference   Contains the transcripts (fasta format)  
 --config      Is the configuration for bowtie (options desired)
 --outdir      Is your output dir
 --snps        If you want to generate a file useful for finding snps (optional)

bowtie must be in the path

The config must specify parameters for bowtie.

=head1 DESCRIPTION

Run and parse bowtie

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

my ($help, $index, $sequence, $config, $outdir, $reference, $snps);
GetOptions(
    'h|help'          => \$help,
    'i|index=s'       => \$index,
    's|sequence=s'    => \$sequence,
    'r|reference=s'   => \$reference,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    'x|snps'        => \$snps,	   
    ) or pod2usage;

pod2usage if $help;

for my $option ($index, $sequence, $outdir, $config, $reference){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/results`;
`mkdir -p $outdir/stderr`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# store the lengths of all transcripts
warn "Storing transcript lengths\n";
my $txdata = {};
my $txin = Bio::SeqIO->new (-format=>'fasta', -file=>"$reference");
while (my $tx_obj = $txin->next_seq()){
    my $id = $tx_obj->display_id();
    my $len = $tx_obj->length();
    $txdata->{$id} = $len;
}


# do the bowtie run
warn "bowtie $sequence\n";
my $alnfile = bowtie_run ($sequence, $index, $conf, $outdir);
#my $alnfile = "s_4_sequence.fastq.alns";
# do the bowtie parse
warn "parse $sequence\n";
bowtie_parse ($alnfile, $txdata, $outdir, $snps);


#####SUBS#####

sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);

    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}


sub bowtie_parse{
    my $alns = shift;
    my $txdata = shift;
    my $out  = shift;
    my $snps = shift;

    # mine the alignment data
    my $mappedrefs = {};
    my $mappedreadcount = 0;
    my $readcache = "tattie";
    open (A, "$out/results/$alns");
    (open (S, ">$out/results/$alns.forsnps")) if ($snps);
    while (my $line = <A>){
	chomp $line;
	
	my ($readname, $strand, $refname, $posref, $readseq, $readquals, $other, $mismatches) = 
	    split (/\t/, $line);
	
	# store and count only the best (first) from the set
	if ($readname ne $readcache){
	    $mappedrefs->{$refname}++;
	    $mappedreadcount++;
	    
	    # for the snp calling file
	    if ($snps){
		if ($strand eq "-"){
		    my $revcompreadseq = revdnacomp ($readseq);
		    my $revcompreadquals = reverse ($readquals);
		    print S "\@$refname\t$readname\n$revcompreadseq\n+$refname\t$readname\n$revcompreadquals\n";
		}
		else {
		    print S "\@$refname\t$readname\n$readseq\n+$refname\t$readname\n$readquals\n";
		}
	    }
	}
	$readcache = $readname;
	
    }
    close (A);
    (close (S)) if ($snps);

    # print out absolute read mappings per reference and rpkm of each reference
    open (R, ">$outdir/results/$alns.parse");
    foreach my $contig (sort keys %$mappedrefs){
	my $rpkm = ($mappedrefs->{$contig} / ($txdata->{$contig} / 1000) / ($mappedreadcount / 1000000));
	print R "$contig\t$mappedrefs->{$contig}\t$rpkm\n";
    }
    close (R);
}

sub bowtie_run{
    my $seq   = shift;
    my $index = shift;
    my $conf  = shift;
    my $out   = shift;

    # get the sequence name
    my $seq_name;
    if ($seq =~/\//g){
        $seq =~m/.*\/(.*)$/;
        $seq_name = $1;
    }

    else {
        $seq_name = $seq;
    }

    # run bowtie
    my $bowtie = "bowtie";
    ($bowtie .= " $conf->{'PARAMETERS'}") if ($conf->{'PARAMETERS'});
    $bowtie .= " $index";
    $bowtie .= " $seq";
    $bowtie .= " 1>$out/results/$seq_name.alns";
    
    `$bowtie 2>$outdir/stderr/$seq_name.stderr`;
    return ("$seq_name.alns");
}

sub parse_config{
    my $file = shift;

    # parse the config file
    my %config;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }
    close (F);

    return (\%config);
}
