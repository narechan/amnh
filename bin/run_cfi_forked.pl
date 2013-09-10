#!/usr/bin/perl -w

=head1 NAME

run_cfi_forked.pl

=head1 SYNOPSIS

  run_cfi_forked.pl

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

PAUP must be in your path
Requires Cfi.pm
Requires the bioperl libraries

The config must specify all parameters for the
search (treecommand), the reps (samples), 
the max trees searched (maxtrees).

Maxtrees can be left undefined in which case it will be
set to 100 (default) and the command, "SET increase=auto" 
will be used.
Examp:

TREECOMMAND=hsearch swap=nni addseq=closest
MAXTREES=200
SAMPLES=10

=head1 DESCRIPTION

Calculate CFI given successive additions of all partitions,
one at a time

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2009 American Museum of Natural History

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
    `cfi_setup.pl -m $matrixfile -c $configfile -o $outdir -s $step -x`;
}
else{
    `cfi_setup.pl -m $matrixfile -c $configfile -o $outdir -s $step`;
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
    warn "CFI run $expt\n";
    `cfi_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$expt`;
    $pm->finish;
}                                                                                
$pm->wait_all_children;    
    
# run the topos script
warn "CFI topologies\n";
`cfi_topos.pl -m $matrixfile -c $configfile -o $outdir`;

# Do the RPAA analysis for one input tree
# or all trees in the analysis
if ($treefile){
    warn "CFI parse $treefile\n";
    `cfi_parse.pl -m $matrixfile -t $treefile -c $configfile -o $outdir`;
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
	warn "CFI parse $index\n";
	`cfi_parse.pl -m $matrixfile -c $configfile -o $outdir -i $index -d $outdir/topology.distribution`;
	$tpm->finish;
    }
    $tpm->wait_all_children;
}
