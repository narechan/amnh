#!/usr/bin/perl -w

=head1 NAME

run_lild_forked.pl

=head1 SYNOPSIS

run_lild_forked.pl

Options:

--procs        is the number of procs designated
--matrix       is the data matrix
--tree         is the total evidence tree    
--config       is your paup config
--outdir       is your outdir for the run
--window       is your window size if doing                                                               
                LILD sildeRule (optional)                                                              
--motion       is the number of residues the                                                         
                window moves per calculation (5' -> 3')                                             
                if doing LILD slideRule (optional)   
--start        is the start position of slideRule (optional)
--end          is the end position of slideRule    (optional)              

Both the matrix and the tree must be in
nexus format.

Requires the bioperl libs.
Requires Lild.pm
Requires that lild_setup.pl, lild_run.pl, and lild_parse.pl
   be in your path.

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

my ($help, $treefile, $matrixfile, $configfile, $procs, $outdir, $window, $motion, $start, $end);
GetOptions(
    'h|help'           => \$help,
    't|tree=s'         => \$treefile,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'p|procs=s'        => \$procs,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $matrixfile, $configfile, $procs, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# run the setup script
warn "Running the setup\n";
if ( ($window) and ($motion) and ($start) and ($end) ){
    `lild_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end`;
}
else {
    `lild_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir`;
}

# run the run script for each experiment
opendir (E, "$outdir/expts");
my @expts = grep (/^.+\..+$/, readdir(E));
closedir (E);

my $pm = Parallel::ForkManager->new($procs);
foreach my $expt (@expts){
    $pm->start and next;
    warn "Running LILD for $expt\n";
    `lild_run.pl -m $matrixfile -c $configfile -o $outdir -t $treefile -e $outdir/expts/$expt`;
    $pm->finish;
}
$pm->wait_all_children;

# run the parse
`lild_parse.pl -m $matrixfile -c $configfile -o $outdir -t $treefile`;
