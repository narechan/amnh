#!/usr/bin/perl -w

=head1 NAME

run_jk_forked.pl

=head1 SYNOPSIS

run_jk_forked.pl

Options:

--procs        is the number of procs designated
--matrix       is the data matrix
--config       is your paup config
--outdir       is your outdir for the run

The matrix must be in
nexus format.

Requires the bioperl libs.
Requires Jk.pm.
Requires paup.
Requires that jk_setup.pl, jk_run.pl, and jk_parse.pl
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

my ($help, $matrixfile, $configfile, $procs, $outdir);
GetOptions(
    'h|help'           => \$help,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'p|procs=s'        => \$procs,
    'o|outdir=s'       => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $configfile, $procs, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# run the setup script
warn "Running the setup\n";
`jk_setup.pl -m $matrixfile -c $configfile -o $outdir`;


# run the run script for each experiment
opendir (E, "$outdir/cmds");
my @expts = grep (/^.+\..+$/, readdir(E));
closedir (E);

my $pm = Parallel::ForkManager->new($procs);
foreach my $expt (@expts){
    $pm->start and next;
    warn "JK for $expt\n";
    `jk_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$expt`;
    $pm->finish;
}
$pm->wait_all_children;

# run the parse
`jk_parse.pl -m $matrixfile -c $configfile -o $outdir`;
