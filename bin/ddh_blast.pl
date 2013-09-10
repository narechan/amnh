#!/usr/bin/perl -w

=head1 NAME

ddh.pl

=head1 SYNOPSIS

ddh.pl

Options:

--procs is the number of || processes
--config is your run configuration
--query is your query genome
--reference is your reference genome
--outdir is your output dir

Requires the bioperl libs. 
Requires Mummer.pm, Metasim.pm, Blast.pm
Requires blast, formatdb, mummer and metasim
    to be in your path

Note that coverages in the config file assume:
454: 280bp; Sanger: 1000bp; solexa: 36bp

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
use Metasim;

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
`mkdir -p $outdir/metasim/reads`;
`mkdir -p $outdir/metasim/stderr`;
`mkdir -p $outdir/blast/raw`;
`mkdir -p $outdir/blast/parse`;
`mkdir -p $outdir/blast/results`;
`mkdir -p $outdir/blast/blastdb`;
`mkdir -p $outdir/summaries`;

#####MAIN#####

# instantiate all objects and load config files
my $mummerobj = Mummer->new;
$mummerobj->load_config ($configfile);
my $blastobj = Blast->new;
$blastobj->load_config ($configfile);
my $metasimobj = Metasim->new;
$metasimobj->load_config ($configfile);

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

# determine the number of reads needed to approx covs
# in the config file and generate the reads
# with metasim in proportion to the size
# of each chrom or plasmid
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

#    last;
    $metasimpm->finish;
}
$metasimpm->wait_all_children;

# get the reads generated by metaim
opendir (R, "$outdir/metasim/reads");
my @queries = grep (/^.+\..+$/, readdir(R));
closedir (R);

# generate blastdb of the ref genome
`cp $reffile $outdir/blast/blastdb`;

my $ref_name;
if ($reffile =~/\//g){
    $reffile =~m/.*\/(.*)$/;
    $ref_name = $1;
}
else {
    $ref_name = $reffile;
}

$blastobj->generate_blastdb("$outdir/blast/blastdb/$ref_name", "$outdir/blast");
`rm $outdir/blast/blastdb/$ref_name`;

# map the metasim generated reads to the ref genome
my $blastpm = Parallel::ForkManager->new($procs);
foreach my $query (sort @queries){
    $blastpm->start and next;

    # analyze the query for true count and len
    my $seqobj = Bio::SeqIO->new(-file   =>"$outdir/metasim/reads/$query",
                                 -format =>"fasta");
    my $counter = 0;
    my $len;
    while (my $seq = $seqobj->next_seq()){
        $counter++;
        $len += $seq->length();
    }

    # start the summary file for this query
    open (SUMF, ">$outdir/summaries/$counter.summary");
    print SUMF "$counter\t";

    my $tcov = $len / $qlen;
    print SUMF "$len\t$tcov\t";

    # print the mummer data
    print SUMF "$mumcov\t$mumcovpct\t";
    
    # run and parse blast
    warn "Blast $queryfile $reffile $query\n";
    my $name = $blastobj->run_blast
	("$outdir/metasim/reads/$query", "$outdir/blast/blastdb/$ref_name", "$outdir/blast/raw");
    warn "Parse $queryfile $reffile $query\n";
    my ($refdata, $qrydata) = $blastobj->parse_blast     # note: leaving identity cutoff at 0
	("$outdir/blast/raw/", $name, "$outdir/blast/parse", 0); 

    # get blast ref coverage
    my $blastcov = $blastobj->genome_coverage ($refdata);
    
    # print blast summary files and 
    # print required data to the overall summary
    my $cover = print_coverage ("$outdir/blast/results/$name.refhits", $blastcov);

    open (SUM, ">$outdir/blast/results/$name.refsum");
    my $refcov;
    foreach my $chrom (sort keys %$cover){
	$refcov += $cover->{$chrom};
	print SUM "$chrom\t$cover->{$chrom}\n";
    }
    close (SUM);
    my $refcovpct = $refcov / $rlen;
    print SUMF "$refcov\t$refcovpct\t";

    open (QRY, ">$outdir/blast/results/$name.qryhits");
    my $qryhits; # measures the number reads that hit somewhere
    my $refhits; # measures the number of ref seqs hit
    my $hsphits; # measures the number of hsps
    foreach my $qry (sort keys %$qrydata){
	$qryhits++;
	foreach my $hit (sort {$a <=> $b} keys %{$qrydata->{$qry}}){
	    $refhits++;
	    $hsphits += $qrydata->{$qry}->{$hit};
	    print QRY "$qry\t$hit\t$qrydata->{$qry}->{$hit}\n";
	}
    }
    close (QRY);
    
    my $fracreads = $qryhits / $counter;
    print SUMF "$qryhits\t$refhits\t$hsphits\t$fracreads\n";
    close (SUMF);

    # purge or compress blast data
    `rm $outdir/blast/raw/$query.out`;
    `gzip $outdir/blast/parse/$query.out.parse`;

    $blastpm->finish
}
$blastpm->wait_all_children;

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

