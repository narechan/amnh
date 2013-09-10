#!/usr/bin/perl -w

=head1 NAME

lild_run.pl

=head1 SYNOPSIS

lild_run.pl

Options:

--experiment is your charset/node pair file
--config is your run configuration
--matrix is your alignfile
--tree is your treefile    
--outdir is your outdir

The matrix must be in nexus format.

Requires the bioperl libs. 
Requires Lild.pm

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
use Lild;

my ($help, $exptfile, $matrixfile, $configfile, $treefile, $outdir);
GetOptions(
    'h|help'          => \$help,
    'e|experiment=s'  => \$exptfile,
    'm|matrix=s'      => \$matrixfile,
    'c|config=s'      => \$configfile,
    't|tree=s'        => \$treefile,   
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($exptfile, $matrixfile, $configfile, $treefile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $lildobj = Lild->new;
$lildobj->load_config ($configfile);
$lildobj->load_aln    ($matrixfile);
$lildobj->load_tree   ($treefile);

# create the paup command file
my ($partition, $node) =
    $lildobj->build_lild (
			  $exptfile, 
			  $matrixfile, 
			  "$outdir/cmds",
			  "$outdir/logs"
			  );

# run paup
$lildobj->run_paup (
		    "$outdir/cmds",
		    $partition,
		    $node
		    );
