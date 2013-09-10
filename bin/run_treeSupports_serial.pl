#!/usr/bin/perl -w

=head1 NAME

run_treeSupports_serial.pl

=head1 SYNOPSIS

run_treeSupports_serial.pl

Options:

--matrix       is the data matrix
--tree         is the total evidence tree    
--config       is your paup config
--outdir       is your outdir

Both the matrix and the tree must be in
nexus format.

Requires the bioperl libs.
Requires TreeSupports.pm
Requires that treeSupports_setup.pl, treeSupports_run.pl, and treeSupports_parse.pl
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
`treeSupports_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir`;

# run the run script for each PBS experiment                                                   
opendir (PBS, "$outdir/pbs/expts");
my @pbsexpts = grep (/^.+\..+$/, readdir(PBS));
closedir (PBS);

foreach my $pbsexpt (@pbsexpts){
    warn "PBS for $pbsexpt\n";
    `treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -s pbs -o $outdir -e $outdir/pbs/expts/$pbsexpt`;
}


# run the run script for each BS experiment
opendir (BS, "$outdir/bs/expts");
my @bsexpts = grep (/^.+\..+$/, readdir(BS));
closedir (BS);

foreach my $bsexpt (@bsexpts){
    warn "BS for $bsexpt\n";
    `treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -s bs -o $outdir -e $outdir/bs/expts/$bsexpt`;
} 

# run the parse
`treeSupports_parse.pl -m $matrixfile -t $treefile -c $configfile -o $outdir`;
