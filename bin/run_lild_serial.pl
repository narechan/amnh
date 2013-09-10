#!/usr/bin/perl -w

=head1 NAME

run_lild_serial.pl

=head1 SYNOPSIS

run_lild_serial.pl

Options:

--matrix       is the data matrix
--tree         is the total evidence tree    
--config       is your paup config
--outdir       is your outdir for the run

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

my ($help, $treefile, $matrixfile, $configfile, $outdir);
GetOptions(
    'h|help'           => \$help,
    't|tree=s'         => \$treefile,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $matrixfile, $configfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# run the setup script
warn "Running the setup\n";
`lild_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir`;

# run the run script for each experiment
opendir (E, "$outdir/expts");
my @expts = grep (/^.+\..+$/, readdir(E));
closedir (E);

foreach my $expt (@expts){
    warn "Running LILD for $expt\n";
    `lild_run.pl -m $matrixfile -c $configfile -o $outdir -t $treefile -e $outdir/expts/$expt`;
} 

# run the parse
`lild_parse.pl -m $matrixfile -c $configfile -o $outdir -t $treefile`;
