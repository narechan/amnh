#!/usr/bin/perl -w

=head1 NAME

insid.pl

=head1 SYNOPSIS

insid.pl is a pipeline script that aligns next generation
query sequence reads at various coverages to either an assembled 
reference genome or an independent set of reference reads. 

Query reads at various coverages can be supplied in the form of 
actual runs of raw nucleotide sequence (note that this program
will not work in colorspace).

The output of this process is a tab-delimited file detailing the
percent of the refernce sequence covered (PRC) by at least one 
read in the query.

This output can then be used to model the saturation kinetics of 
the iterative alignment procedure.  

Dependencies:

Requires the bioperl libraries.
Requires the installation of mummer into $PATH.
Requires the runner module Mummer.pm. 

Options:

--contigs is set to true if input is contigs and not reads
--procs is the number of parallel processes to fork
--config is your run configuration
--reference is your reference genome
--seqin is the dir containing your input query reads or your
    fasta file containing assembled contigs
--outdir is your output dir

Command Line Details:

The config file specifies options for mummer

MUMMER=-mum -b -c -F -l 19

All sequence files should be in fasta format.

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

my ($help, $procs, $configfile, $reffile, $outdir, $seqin, $contigs);
GetOptions(
    'h|help'          => \$help,
    'p|procs=s'       => \$procs,
    'c|config=s'      => \$configfile,
    'r|reference=s'   => \$reffile,	   
    'o|outdir=s'      => \$outdir,
    's|seqin=s'       => \$seqin,
    't|contigs'       => \$contigs,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $procs, $reffile, $outdir, $seqin){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

# make the dirs to hold the data
`mkdir -p $outdir/reads/raw`;
`mkdir -p $outdir/reads/results`;
`mkdir -p $outdir/reads/stderr`;
`mkdir -p $outdir/summaries`;

#####MAIN#####

# instantiate all objects and load config files
my $mummerobj = Mummer->new;
$mummerobj->load_config ($configfile);

# date and time trackers
my $date;
my $time;

# get the query provided at command line
($date, $time) = time_stamp();
warn "PRC start $date $time\n";

my @queries;
if ($contigs){
    push (@queries, $seqin);
    $procs = 1; #force one fork
}
else{
    opendir (R, "$seqin");
    @queries = grep (/^.+\..+$/, readdir(R));
    closedir (R);
}

# count the number of queries
my $queries = @queries;

my $rlen;
my $rseqobj = Bio::SeqIO->new(-file   =>"$reffile",
			      -format =>"fasta");
while (my $rseq = $rseqobj->next_seq()){
    $rlen += $rseq->length;
}

# map the reads to the ref genome or reference reads
my $mummerpm = Parallel::ForkManager->new($procs);
foreach my $query (sort @queries){
    $mummerpm->start and next;
    
    # analyze the query for true count and len                                    
    my $seqobj;
    if ($queries > 1){
	$seqobj = Bio::SeqIO->new(-file   =>"$seqin/$query",
				  -format =>"fasta");
    }
    else {
	$seqobj = Bio::SeqIO->new(-file   =>"$query",
				  -format =>"fasta");
    }
    
    my $counter = 0;
    my $qlen;
    while (my $seq = $seqobj->next_seq()){
        $counter++;
        $qlen += $seq->length();
    }
    
    # run and parse reads mummer                                                                           
    `mkdir -p $outdir/reads/raw/$query`;
    `mkdir -p $outdir/reads/results/$query`;
    `mkdir -p $outdir/reads/stderr/$query`;
    
    if ($queries > 1){
	$mummerobj->run_mummer
	    ("$seqin/$query", $reffile, "$outdir/reads/raw/$query", "$outdir/reads/stderr/$query");
    }
    else{
	$mummerobj->run_mummer
            ($query, $reffile, "$outdir/reads/raw/$query", "$outdir/reads/stderr/$query");
    }

    my $mummerdata = $mummerobj->parse_mummer ("$outdir/reads/raw/$query/out.mummer");

    # if there were alignments, analyze them
    my $refcov;
    if (@$mummerdata){
	
	# determine reads mummer ref coverage                                                             
	my $mummercov = $mummerobj->genome_coverage ($mummerdata);
	
	# print out the reads mummer summary files                                                       
	my $covered = print_coverage ("$outdir/reads/results/$query/out.refhits", $mummercov);
	
	open (SUM, ">$outdir/reads/results/$query/out.refsum");
	foreach my $chrom (sort keys %$covered){
	    $refcov += $covered->{$chrom};
	    print SUM "$chrom\t$covered->{$chrom}\n";
	}
	close (SUM);
    }
    else {
	$refcov = 1; #small but non-zero, otherwise insidR.pl (mm.R) will break
    }
    
    # start the summary file for this query
    my $refcovpct = $refcov / $rlen;
    open (SUMF, ">$outdir/summaries/$counter.summary");
    print SUMF "$counter\t$qlen\t";
    print SUMF "$rlen\t$refcov\t$refcovpct\n";
    close (SUMF);

    ($date, $time) = time_stamp();
    warn "Mummer end $query $reffile $date $time\n";
    
    $mummerpm->finish
    }
$mummerpm->wait_all_children;

`cat $outdir/summaries/* | sort -n > $outdir/summary`;

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

sub time_stamp {
    my ($d,$t);
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

    $year += 1900;
    $mon++;
    $d = sprintf("%4d-%2.2d-%2.2d",$year,$mon,$mday);
    $t = sprintf("%2.2d:%2.2d:%2.2d",$hour,$min,$sec);
    return($d,$t);
}
