#!/usr/bin/perl -w

=head1 NAME

jk_run.pl

=head1 SYNOPSIS

jk_run.pl

Options:

--expt is the expt / cmd you want to run
--config is your run configuration
--matrix is your alignfile
--outdir is your outdir

The matrix must be in nexus format.

Requires the bioperl libs. 
Requires Jk.pm
Requires paup

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
use Jk;

my ($help, $exptfile, $matrixfile, $configfile, $outdir);
GetOptions(
    'h|help'          => \$help,
    'e|experiment=s'  => \$exptfile,
    'm|matrix=s'      => \$matrixfile,
    'c|config=s'      => \$configfile,
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($exptfile, $matrixfile, $configfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $jkobj = Jk->new;
$jkobj->load_config ($configfile);

# run paup
$jkobj->run_paup ($exptfile);
