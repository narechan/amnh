#!/usr/bin/perl -w

=head1 NAME

insid.pl

=head1 SYNOPSIS

insid.pl is a pipeline script that aligns next generation
query sequence reads at various coverages to either an assembled 
reference genome or an independent set of reference reads. 

Query reads at various coverages can be supplied in the form of 
actual runs of illumina or 454 sequence, or can be simulated from an
assembled query genome.

The output of this process is a tab-delimited file detailing the
percent of the refernce sequence covered (PRC) by at least one 
reads in the query.

This output can then be used to model the saturation kinetics of 
the iterative alignment procedure.  

Dependencies:

Requires the bioperl libraries.
Requires the installation of mummer and metasim into $PATH.
Requires the runner modules Mummer.pm and Metasim.pm. 

Options:

--procs/-p is the number of parallel processes to fork
--config/-c is your run configuration
--reference/-r is your reference genome
--seqin/-s is the dir containing your input query sequences (optional)
    the sequences option will short circuit metasim and 
    take reads from the input dir for analysis
--query/-q is your query genome (optional)
--outdir/-o is your output dir

Command Line Details:

The config file specifies options for mummer and metasim, and coverages
for query read simulation (if applicable). As in:

MUMMER=-mum -b -c -F -l 19
METASIM=--empirical
COVERAGES=0.05,0.1,0.25,0.5,0.75,1,...

For the query, either a query genome or seqin directory must be
specified.

All sequence files should be in fasta format.

Coverages in the config file assume the following sequence lengths:
454: 280bp; Sanger: 1000bp; solexa: 36bp

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2010 American Museum of Natural History

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
use Metasim;

my ($help, $procs, $configfile, $queryfile, $reffile, $outdir, $seqin);
GetOptions(
    'h|help'          => \$help,
    'p|procs=s'       => \$procs,
    'c|config=s'      => \$configfile,
    'q|query=s'       => \$queryfile,
    'r|reference=s'   => \$reffile,	   
    'o|outdir=s'      => \$outdir,
    's|seqin=s'       => \$seqin,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $procs, $reffile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

# make the dirs to hold the data
`mkdir -p $outdir/genome/raw`;
`mkdir -p $outdir/genome/results`;
`mkdir -p $outdir/genome/stderr`;
`mkdir -p $outdir/metasim/reads`;
`mkdir -p $outdir/metasim/stderr`;
`mkdir -p $outdir/reads/raw`;
`mkdir -p $outdir/reads/results`;
`mkdir -p $outdir/reads/stderr`;
`mkdir -p $outdir/summaries`;

#####MAIN#####

# instantiate all objects and load config files
my $mummerobj = Mummer->new;
$mummerobj->load_config ($configfile);

my $metasimobj = Metasim->new;
$metasimobj->load_config ($configfile);

# find the total length of the query and ref seqs                                                         
# unless already specifying multiple read set queries
my $qlen;
my $rlen;
my $mumcov;
my $mumcovpct;
unless ($seqin){
    my $qseqobj = Bio::SeqIO->new(-file   =>"$queryfile",
				  -format =>"fasta");
    while (my $qseq = $qseqobj->next_seq()){
	$qlen += $qseq->length;
    }
    
    my $rseqobj = Bio::SeqIO->new(-file   =>"$reffile",
				  -format =>"fasta");
    while (my $rseq = $rseqobj->next_seq()){
	$rlen += $rseq->length;
    }
    
    # run and parse genome mummer
    warn "Mummer $queryfile $reffile\n";
    $mummerobj->run_mummer
	($queryfile, $reffile, "$outdir/genome/raw", "$outdir/genome/stderr");
    my $mummerdata = $mummerobj->parse_mummer ("$outdir/genome/raw/out.mummer");
    
    # determine genome mummer ref coverage
    my $mummercov = $mummerobj->genome_coverage ($mummerdata);
    
    # print out the genome mummer summary files
    my $covered = print_coverage ("$outdir/genome/results/out.refhits", $mummercov);
    
    open (SUM, ">$outdir/genome/results/out.refsum");
    foreach my $chrom (sort keys %$covered){
	$mumcov += $covered->{$chrom};
	print SUM "$chrom\t$covered->{$chrom}\n";
    }
    close (SUM);
    $mumcovpct = $mumcov / $rlen;
    
    # determine the number of reads needed to approx covs
    # in the config file and generate the reads
    # with metasim in proportion to the size
    # of each chrom or plasmid.
    # short circuit of sequences already provided
    my $metasimpm = Parallel::ForkManager->new($procs);
    foreach my $cov (@{$metasimobj->get_coverages}){
	my $covname = $cov; 
	$covname =~s/\./-/;
	$metasimpm->start and next;
	
	my $type = $metasimobj->get_metasim_params;
	my $readlen;
	if ($type =~m/454/){
	    $readlen = 280;
	}
	elsif ($type =~m/sanger/){
	    $readlen = 1000;
	}
	elsif ($type =~m/empirical/){
	    $readlen = 36;
	}
	else {
	    print STDERR "Unknown readlen for seqtype\n";
	    die;
	}
	
	my $reads = int( (($cov * $qlen) / $readlen) + 0.5 );
	($reads = int (($reads / 2) + 0.5)) if ($type =~m/sanger/); # HACK! compensates for bug in metasim
	
	warn "Metasim $queryfile $reffile $cov $readlen\n";
	`mkdir -p $outdir/metasim/reads/$covname`;
	`mkdir -p $outdir/metasim/stderr/$covname`;
	my $seqobj = Bio::SeqIO->new(-file   =>"$queryfile",
				     -format =>"fasta");
	my $counter = 0;
	while (my $seq = $seqobj->next_seq()){
	    $counter++;
	    my $id       = $seq->display_id();
	    my $sequence = $seq->seq();
	    my $len      = $seq->length();
	    my $seqreads = int ( ($reads * ($len / $qlen)) + 0.5 );
	    
	    # skip if we don't want any reads from the seq
	    next if ($seqreads == 0);
	    
	    # skip if the clone size is too small
	    # using 5x the avg clone size (2000 for 454, 5000 for sanger)
	    next if ( ($len < 10000) and ($type =~m/454/) );
	    next if ( ($len < 25000) and ($type =~m/sanger/) );
	    
	    open (TC, ">$outdir/metasim/reads/$covname/temp.fasta");
	    print TC ">$id\n$sequence\n";
	    close (TC);
	    
	    $metasimobj->generate_reads 
		($seqreads, 
		 "$outdir/metasim/reads/$covname/temp.fasta", 
		 "$outdir/metasim/reads/$covname", 
		 "$outdir/metasim/stderr/$covname",
		 $counter);
	    
	    `rm $outdir/metasim/reads/$covname/temp.fasta`;
	}
	`cat $outdir/metasim/reads/$covname/*.fasta > $outdir/metasim/reads/$reads.fasta`;
	`cat $outdir/metasim/stderr/$covname/stderr.* > $outdir/metasim/stderr/$reads.stderr`;
	`rm -rf $outdir/metasim/reads/$covname`;
	`rm -rf $outdir/metasim/stderr/$covname`;
	
	$metasimpm->finish;
    }
    $metasimpm->wait_all_children;
}

# get the reads generated by metasim or provided in command line
my @queries;
my $querystem;
if ($seqin){
    opendir (R, "$seqin");
    @queries = grep (/^.+\..+$/, readdir(R));
    $querystem = $seqin;
    closedir (R);
    
    my $rseqobj = Bio::SeqIO->new(-file   =>"$reffile",
				  -format =>"fasta");
    while (my $rseq = $rseqobj->next_seq()){
        $rlen += $rseq->length;
    }
}
else{
    opendir (R, "$outdir/metasim/reads");
    @queries = grep (/^.+\..+$/, readdir(R));
    $querystem = "$outdir/metasim/reads";
    closedir (R);
}

# map the reads to the ref genome or reference reads
my $mummerpm = Parallel::ForkManager->new($procs);
foreach my $query (sort @queries){
    $mummerpm->start and next;
    
    # analyze the query for true count and len                                    
    my $seqobj = Bio::SeqIO->new(-file   =>"$querystem/$query",
                                 -format =>"fasta");
    my $counter = 0;
    my $len;
    while (my $seq = $seqobj->next_seq()){
        $counter++;
        $len += $seq->length();
    }
    
    # run and parse reads mummer                                                                           
    warn "Mummer $query $reffile\n";
    `mkdir -p $outdir/reads/raw/$query`;
    `mkdir -p $outdir/reads/results/$query`;
    `mkdir -p $outdir/reads/stderr/$query`;
    
    $mummerobj->run_mummer
	("$querystem/$query", $reffile, "$outdir/reads/raw/$query", "$outdir/reads/stderr/$query");
    my $mummerdata = $mummerobj->parse_mummer ("$outdir/reads/raw/$query/out.mummer");
    
    # determine reads mummer ref coverage                                                                     
    my $mummercov = $mummerobj->genome_coverage ($mummerdata);

    # print out the reads mummer summary files                                                       
    my $covered = print_coverage ("$outdir/reads/results/$query/out.refhits", $mummercov);
    
    my $refcov;
    open (SUM, ">$outdir/reads/results/$query/out.refsum");
    foreach my $chrom (sort keys %$covered){
	$refcov += $covered->{$chrom};
	print SUM "$chrom\t$covered->{$chrom}\n";
    }
    close (SUM);
    
   # start the summary file for this query
    my $tcov;
    my $refcovpct = $refcov / $rlen;
    open (SUMF, ">$outdir/summaries/$counter.summary");
    unless ($seqin){
	$tcov = $len / $qlen;
	print SUMF "$counter\t$len\t$tcov\t";
	print SUMF "$mumcov\t$mumcovpct\t";
	print SUMF "$rlen\t$refcov\t$refcovpct\n";
    }
    else {
	print SUMF "$counter\t$len\tUNKNOWN\t"; #query coverage unknown b/c query genome size unknown
	print SUMF "UNKNOWN\tUNKNOWN\t"; #no genome level mummer data available when query genome not supplied
	print SUMF "$rlen\t$refcov\t$refcovpct\n";
    }
    close (SUMF);
    
    $mummerpm->finish
    }
$mummerpm->wait_all_children;

`tar czf $outdir/metasim.tar.gz $outdir/metasim`;
`tar czf $outdir/reads.tar.gz $outdir/reads`;
`rm -rf $outdir/metasim`;
`rm -rf $outdir/reads`;
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

