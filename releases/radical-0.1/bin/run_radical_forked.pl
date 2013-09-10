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
 --procs       Is the number of procs for fork (if forking)
 --step        Is the stepsize across partitions (i.e., 5 to run every 5th partition)                  
                 (optional)
 --subx        Is set only if you want sub-matrices for each section of the curve
  
The aln file must be in nexus format.

Dependencies:

PAUP or RAxML must be installed and in $PATH.
Requires Radical.pm
Requires the bioperl libraries
Requires Statistics::Descriptive
Requires Parallel:ForkManager

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

my ($help, $matrixfile, $treefile, $outdir, $configfile, $procs, $step, $subx);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    't|tree=s'        => \$treefile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    'p|procs=s'       => \$procs,
    's|step=s'        => \$step,
    'x|subx'          => \$subx,
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

# run the run script with each expt
opendir (E, "$outdir/cmds");
my @expts = sort (readdir(E));
shift @expts;
shift @expts;
closedir (E);

my $pm = Parallel::ForkManager->new($procs);                                                    
foreach my $expt (@expts){
    $pm->start and next;
    warn "RADICAL run $expt\n";
    `radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$expt`;
    $pm->finish;
}                                                                                
$pm->wait_all_children;    
    
# run the topos script
warn "RADICAL topologies\n";
`radical_topos.pl -m $matrixfile -c $configfile -o $outdir`;

# Do the RADICAL analysis for one input tree
# or all trees in the analysis
if ($treefile){
    warn "RADICAL parse $treefile\n";
    `radical_parse.pl -m $matrixfile -t $treefile -c $configfile -o $outdir`;
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
	warn "RADICAL parse $index\n";
	`radical_parse.pl -m $matrixfile -c $configfile -o $outdir -i $index -d $outdir/topology.distribution`;
	$tpm->finish;
    }
    $tpm->wait_all_children;
}
