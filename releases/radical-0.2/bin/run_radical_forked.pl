#!/usr/bin/perl -w

=head1 NAME

run_radical_forked.pl

=head1 SYNOPSIS

  run_radical_forked.pl -- Run the RADICAL analysis on multicore architecture.

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --tree        Is your TE treefile (optional)
 --outdir      Is your output dir
 --config      Is the configuration for your tests (sections can be a single number for divisor, or
    a comma sep list of concat points desired)
 --procs       Is the number of procs for fork
 --step        Is the stepsize across partitions (i.e., 5 to run every 5th partition) (optional)
 --subx        Is set only if you want sub-matrices for each section of the curve (optional)
 --checking    Is set if you want to abort concatenation when the statistical criterion has been
    reached for the replicates. If you want to use the abort procedure, a treefile must be specified (optional)

The aln file must be in nexus format.

Dependencies:

PAUP or RAxML must be installed and in $PATH.
Requires Radical.pm
Requires the bioperl libraries
Requires Statistics::Descriptive
Requires Parallel:ForkManager
Requires Bio::Phylo

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2011 American Museum of Natural History

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

my ($help, $matrixfile, $treefile, $outdir, $configfile, $procs, $step, $subx, $checking);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    't|tree=s'        => \$treefile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    'p|procs=s'       => \$procs,
    's|step=s'        => \$step,
    'x|subx'          => \$subx,
    'z|checking'      => \$checking, 
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrixfile, $outdir, $configfile, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# run the setup script
# put step to default 1 if not specified                                    
($step = 1) unless ($step);
warn "CFI setup\n";
if ($subx){
    `radical_setup.pl -m $matrixfile -c $configfile -o $outdir -s $step -x`;
}
else{
    `radical_setup.pl -m $matrixfile -c $configfile -o $outdir -s $step`;
}

# get all the experiments
opendir (E, "$outdir/cmds");
my @exptfiles = sort (readdir(E));
shift @exptfiles;
shift @exptfiles;
closedir (E);

my %samps;
my %expts;
foreach my $file (@exptfiles){
    my ($samp, $expt, $nex) = split (/\./, $file);
    $samps{$samp} = 1;
    $expts{$expt} = 1;
}

my $samps = keys %samps;
my $expts = keys %expts;
my @expts = sort {$a <=> $b} keys %expts; # using array in case step defined     

# execute jobs with checking or do the entire set
if ($checking){
    foreach my $expt (@expts){
	my $pm = Parallel::ForkManager->new($procs);
	foreach my $samp (sort {$a <=> $b} keys %samps){
	    $pm->start and next;
	    warn "RADICAL run $samp.$expt.cmds\n";
	    `radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$samp.$expt.cmds`;
	    $pm->finish;
	}
	$pm->wait_all_children;
	
	warn "Checking Stats\n";
	`radical_check.pl -m $matrixfile -t $treefile -c $configfile -o $outdir > $outdir/temp.stats`;
	
	open (S, "$outdir/temp.stats");
        my $state = <S>;
        chomp $state;
        warn "$state\n";

	if ($state =~m/^STABLE/){
            `rm $outdir/temp.stats`;
            last;
        }
        else {
            `rm $outdir/temp.stats`;
            next;
        }
    }
}
else{
    my $pm = Parallel::ForkManager->new($procs);                                                    
    foreach my $expt (@exptfiles){
	$pm->start and next;
	warn "RADICAL run $expt\n";
	`radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$expt`;
	$pm->finish;
    }                                                                                
    $pm->wait_all_children;    
}

# run the topos script
warn "RADICAL topologies\n";
opendir (T, "$outdir/trees");
my @treefiles = sort (readdir(T));
shift @treefiles;
shift @treefiles;
closedir (T);

foreach my $tree (@treefiles){
    print STDERR "$tree\n";
    `radical_topos.pl -l $outdir/summaries/taxa -t $outdir/trees/$tree >> $outdir/temp.trees`;
}

my $topologies = {};
open (TF, "$outdir/temp.trees");
while (my $line = <TF>){
    chomp $line;
    $line =~s/\;//;
    $topologies->{$line}++;
}
close (TF);

open (TS, ">$outdir/topology.distribution");
my $tcounter = 0;
foreach my $topo (keys %$topologies){
    $tcounter++;
    print TS "$tcounter\t$topo\t$topologies->{$topo}\n";
}
close (TS);

`rm $outdir/temp.trees`;

# Do the RADICAL post process analysis for one input tree
# or all uniq trees generated
if ($treefile){
    warn "RADICAL post-process $treefile\n";
    `radical_pp.pl -m $matrixfile -t $treefile -c $configfile -o $outdir`;
}
else {
    my %topos;
    open (TD, "$outdir/topology.distribution");
    while (my $line = <TD>){
        chomp $line;
	my ($num, $topo) = split (/\t/, $line);
	$topos{$num} = $topo;
    }
    close (TD);

    my $tpm = Parallel::ForkManager->new($procs);
    foreach my $index (keys %topos){
	$tpm->start and next;
	warn "RADICAL post-process $index\n";
	sleep 5; #avoids race to create the consensus trees
	`radical_pp.pl -m $matrixfile -c $configfile -o $outdir -i $index -d $outdir/topology.distribution`;
	$tpm->finish;
    }
    $tpm->wait_all_children;
}

