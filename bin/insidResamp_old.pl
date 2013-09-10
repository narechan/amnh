#!/usr/bin/perl -w

=head1 NAME

insidResamp.pl

=head1 SYNOPSIS

The script wraps around insidP.pl, insid.pl, insidR.pl, and insidD.pl
to either bootstrap or jackknife an insid reads tree. It iterates over the 
entire insid pipeline for every replicate and generates a majority 
rules bootstrap or jackknife tree for the entire process.

Note that for resamp, query coverages must be provided to as total base
counts.

Dependencies:                                                                                            

All insid components and mummer must be in your path
Requires the bioperl libraries.                                                                               
Requires the runner module Mummer.pm.
Requires R.

Options:

--procs is the number of processors to use (note that this parallelizes 
    over the pairwise replicates, not within them)
--data is a directory containing read fasta files for
    all species
--mode is the sampling mode (bootstrap = 1; jackknife = 2)
--config is the mummer configuration file
--refbases is a file containing the base counts                                                            
    for the ref set
--qrybases is a file containing the base counts 
    for the qry bases
--mmr is the path to the mm.R script
--outdir is the output dir for all data
--reps is the number of replicates you want
--clean if you want to remove heavy interim files

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

my ($help, $datadir, $outdir, $refbasefile, $qrybasefile, $mode, $configfile, $mmrscript, $reps, $procs, $clean);
GetOptions(
    'h|help'          => \$help,
    'd|datadir=s'     => \$datadir,
    'o|outdir=s'      => \$outdir,
    'a|refbases=s'    => \$refbasefile,
    'b|qrybases=s'    => \$qrybasefile,
    'c|config=s'      => \$configfile,
    'm|mode=s'        => \$mode,
    's|mmr=s'         => \$mmrscript,
    'r|reps=s'        => \$reps,
    'p|procs=s'       => \$procs,
    'x|clean'       => \$clean,
    ) or pod2usage;
pod2usage if $help;

for my $option ($datadir, $outdir, $refbasefile, $qrybasefile, $mode, $configfile, $mmrscript, $reps, $procs){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir/reps`;

#####MAIN#####    

# get all orgs
opendir (D, "$datadir");
my @readfiles = sort (readdir (D));
shift @readfiles;
shift @readfiles;
closedir (D);

# get the basecounts
my @refbasecounts;
open (B, "$refbasefile");
while (my $basecount = <B>){
    chomp $basecount;
    push (@refbasecounts, $basecount);
}
close (B);

my @qrybasecounts;
open (Q, "$qrybasefile");
while (my $basecount = <Q>){
    chomp $basecount;
    push (@qrybasecounts, $basecount);
}
close (Q);

for (my $i = 1; $i <= $reps; $i++){
    `mkdir -p $outdir/reps/$i/results`;
    `mkdir -p $outdir/reps/$i/models`;
    `mkdir -p $outdir/reps/$i/matrices`;
    
    # resample to create query and ref reads sets for this rep
    my $samplepm = Parallel::ForkManager->new($procs);
    foreach my $readfile (@readfiles){
	$samplepm->start and next;

	print STDERR "INSIDP for $readfile queries\n";
	`insidP.pl -i $datadir/$readfile -o $outdir/reps/$i/qry_reads/$readfile/ -b $qrybasefile -m $mode`;

	print STDERR "INSIDP for $readfile references\n";
	`insidP.pl -i $datadir/$readfile -o $outdir/reps/$i/ref_reads/$readfile/ -b $refbasefile -m $mode`;

	$samplepm->finish;
    }
    $samplepm->wait_all_children;

    # do the insid analysis for every set of query reads against every 
    # reference basecount across all pairwise comparisons in this rep.
    # also do plots and models for each run
    foreach my $basecountR (@refbasecounts){
	foreach my $reffile (@readfiles){

	    my $sampleim = Parallel::ForkManager->new($procs);
	    foreach my $qryfile (@readfiles){
		$sampleim->start and next;
		
		# insid
		print STDERR "INSID for $basecountR $reffile $qryfile\n"; 
		`insid.pl -p 1 -c $configfile -s $outdir/reps/$i/qry_reads/$qryfile -r $outdir/reps/$i/ref_reads/$reffile/$basecountR.fa -o $outdir/reps/$i/results/$basecountR/$qryfile-$reffile`;
		`mv $outdir/reps/$i/results/$basecountR/$qryfile-$reffile/summary $outdir/reps/$i/results/$basecountR/$qryfile-$reffile/$qryfile-$reffile.summary`;
		
		# model
		foreach my $basecountQ (@qrybasecounts){
		    print STDERR "MODEL for $basecountQ $basecountR $reffile $qryfile\n";
		    `insidR.pl -c $basecountQ -r $mmrscript -i $outdir/reps/$i/results/$basecountR/$qryfile-$reffile/$qryfile-$reffile.summary -o $outdir/reps/$i/results/$basecountR/$qryfile-$reffile/model-$basecountQ > $outdir/reps/$i/results/$basecountR/$qryfile-$reffile/model-$basecountQ.out`;
		}
		$sampleim->finish;
	    }
	    $sampleim->wait_all_children;
	} 
	
	# build the matrix for each reference basecount
	foreach my $basecountQ (@qrybasecounts){
	    `find $outdir/reps/$i/results/$basecountR/ -type f -name model-$basecountQ.out -exec cat {} \\; > $outdir/reps/$i/models/$basecountQ-$basecountR.models\n`;
	    `insidD.pl -i $outdir/reps/$i/models/$basecountQ-$basecountR.models -o $outdir/reps/$i/matrices/$basecountQ-$basecountR`;
	}
    }

    # clean up if desired
    if ($clean){
	`rm -rf $outdir/reps/$i/qry_reads/*`;
	`rm -rf $outdir/reps/$i/ref_reads/*`;
	`rm -rf $outdir/reps/$i/results/*`;
    }
}
