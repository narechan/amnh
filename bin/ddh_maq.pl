#!/usr/bin/perl -w

=head1 NAME

ddh_maq.pl

=head1 SYNOPSIS

ddh_maq.pl

Options:

--procs is the number of || processes
--config is your run configuration
--query is your query genome
--reference is your reference genome
--outdir is your output dir

Requires the bioperl libs. 
Requires Mummer.pm, Maq.pm, Blast.pm (for mummer ref covered)
Requires mummer and maq
    to be in your path

Note that coverages in the config file assume:
solexa: 36bp

=head1 DESCRIPTION


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

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

use Bio::SeqIO;

use Mummer;
use Blast;
use Maq;

my ($help, $procs, $configfile, $queryfile, $reffile, $outdir);
GetOptions(
    'h|help'          => \$help,
    'p|procs=s'       => \$procs,
    'c|config=s'      => \$configfile,
    'q|query=s'       => \$queryfile,
    'r|reference=s'   => \$reffile,	   
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $procs, $queryfile, $reffile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

# make the dirs to hold the data
`mkdir -p $outdir/mummer/raw`;
`mkdir -p $outdir/mummer/results`;
`mkdir -p $outdir/mummer/stderr`;
`mkdir -p $outdir/maq/stderr`;
`mkdir -p $outdir/maq/results`;
`mkdir -p $outdir/maq/data`;
`mkdir -p $outdir/summaries`;

#####MAIN#####

# instantiate all objects and load config files
my $mummerobj = Mummer->new;
$mummerobj->load_config ($configfile);
my $blastobj = Blast->new;
$blastobj->load_config ($configfile);
my $maqobj = Maq->new;
$maqobj->load_config ($configfile);

# find the total length of the query and ref seqs                                                         
my $qlen;
my $qseqobj = Bio::SeqIO->new(-file   =>"$queryfile",
                              -format =>"fasta");
while (my $qseq = $qseqobj->next_seq()){
    $qlen += $qseq->length;
}

my $rlen;
my $rseqobj = Bio::SeqIO->new(-file   =>"$reffile",
                              -format =>"fasta");
while (my $rseq = $rseqobj->next_seq()){
    $rlen += $rseq->length;
}

# run and parse mummer
warn "Mummer $queryfile $reffile\n";
$mummerobj->run_nucmer
    ($queryfile, $reffile, "$outdir/mummer/raw", "$outdir/mummer/stderr");
my $mummerdata = $mummerobj->parse_nucmer ("$outdir/mummer/raw/out.showcoords");

# determine mummer ref coverage
my $mummercov = $blastobj->genome_coverage ($mummerdata);

# print out the mummer summary files
my $covered = print_coverage ("$outdir/mummer/results/out.refhits", $mummercov);

my $mumcov;
open (SUM, ">$outdir/mummer/results/out.refsum");
foreach my $chrom (sort keys %$covered){
    $mumcov += $covered->{$chrom};
    print SUM "$chrom\t$covered->{$chrom}\n";
}
close (SUM);
my $mumcovpct = $mumcov / $rlen;

# transform reference into binary                                                                            
my $refname = $maqobj->maq_fasta2bfa
    ($reffile,
     "$outdir/maq/data",
     "$outdir/maq/stderr");

# do a maq expt for every coverage
my $simpm = Parallel::ForkManager->new($procs);
foreach my $cov (@{$maqobj->get_coverages}){
    $simpm->start and next;

    # determine the number of reads needed to approx covs 
    my $readlen  = $maqobj->get_length;
    my $datafile = $maqobj->get_datafile;
    my $reads    = int( (($cov * $qlen) / $readlen) + 0.5 );
    my $pereads  = int ( ($reads / 2) + 0.5 );
    
    # generate the reads
    warn "Maq simulate $queryfile $datafile $cov\n";
    $maqobj->maq_simulate 
	($pereads, 
	 $queryfile, 
	 $datafile,
	 "$outdir/maq/data",
	 "$outdir/maq/stderr");
    
    # transform reads into binary
    for my $pair (1, 2){
	$maqobj->maq_fastq2bfq
	    ("$outdir/maq/data/$pereads-$pair.fastq",
	     "$outdir/maq/data/$pereads-$pair.bfq",
	     "$outdir/maq/stderr/$pereads-$pair.bfq.stderr");
    }
    
    # do the maq match
    warn "Maq match $queryfile $datafile $cov\n";
    $maqobj->maq_match
	($pereads,
	 $refname,
	 "$outdir/maq/data",
	 "$outdir/maq/results",
	 "$outdir/maq/stderr");

    # generate the readable aln data
    warn "Maq process $queryfile $datafile $cov\n";
    $maqobj->maq_mapview
	($pereads,
	 $refname,
	 "$outdir/maq/results",
         "$outdir/maq/stderr");

    # pileup the data per ref position
    $maqobj->maq_pileup
	($pereads,
	 $refname,
	 "$outdir/maq/data",
	 "$outdir/maq/results",
         "$outdir/maq/stderr");

    # parse the pileup file to determine
    # amnt of ref covered
    open (P, "$outdir/maq/results/$pereads-v-$refname.pileup");
    my $maqcov = 0;
    while (my $line = <P>){
	chomp $line;
	my ($acc, $pos, $base, $covered, $data) = 
	    split (/\t/, $line, 5);
	($maqcov++) if ($covered > 0);
    }

    # get the true coverage of the sim reads
    my $simlen = 0;
    for my $pair (1, 2){
	my $seqobj = Bio::SeqIO->new(-file   =>"$outdir/maq/data/$pereads-$pair.fastq",
				     -format =>"fastq");
	while (my $seq = $seqobj->next_seq()){
	    $simlen += $seq->length;
	}
    }

    my $coverage = $simlen / $qlen;
    my $maqcovpct = $maqcov / $rlen;

    open (S, ">$outdir/summaries/$pereads.summary");
    print S "$pereads\t$simlen\t$coverage\t";
    print S "$mumcov\t$mumcovpct\t$maqcov\t$maqcovpct\n";
    close (S);
#    last;
    $simpm->finish;
}
$simpm->wait_all_children;
`cat $outdir/summaries/* > $outdir/summary`;


####SUBS####

sub print_coverage{
    my $file = shift;
    my $covs = shift;

    my $covered ={};
    open (REF, ">$file");
    foreach my $chrom (sort keys %$covs){
	foreach my $region (sort {$a <=> $b} keys %{$covs->{$chrom}}){
	    my $start = $covs->{$chrom}->{$region}[0];
	    my $end   = $covs->{$chrom}->{$region}[1];
	    my $len   = $end - $start + 1;

	    $covered->{$chrom} += $len;
	    print REF "$chrom\t$region\t$start\t$end\n";
	}
    }
    close (REF);

    return ($covered);
}
