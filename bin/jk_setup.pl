#!/usr/bin/perl -w

=head1 NAME

jk_setup.pl

=head1 SYNOPSIS

jk_setup.pl

Options:

--config       is the config file
--matrix       is the data matrix
--outdir       is your output dir for the run

the matrix must be in
nexus format.

Requires the bioperl libs.
Requires Jk.pm

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

my ($help, $configfile, $matrixfile, $outdir);
GetOptions(
    'h|help'           => \$help,
    'c|config=s'       => \$configfile,   
    'm|matrix=s'       => \$matrixfile,
    'o|outdir=s'       => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $matrixfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}


# create dir structure for the results
`mkdir -p $outdir/cmds`;
`mkdir -p $outdir/logs`;

#####MAIN#####

# instantiate the object and load required data
my $jkobj = Jk->new;
$jkobj->load_config ($configfile);

# get the pctdelete's
my $pctdel = $jkobj->get_pctdel;
my @pctdel = split (/\,/, $pctdel);

# generate the paup jk cmds
foreach my $pdel (@pctdel){
    $jkobj->generate_jk_part ($pdel,
			      "$outdir/cmds",
			      "$outdir/logs",
			      $matrixfile);
}

