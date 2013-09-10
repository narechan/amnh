#!/usr/bin/perl -w

=head1 NAME

cfi_run.pl

=head1 SYNOPSIS

  cfi_run.pl

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --outdir      Is your output dir
 --config      Is the configuration for your tests
 --exptfile    Is the expt / cmd you want to run

The aln file must be in nexus format.

PAUP must be in your path
Requires Cfi.pm

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
use Cfi;

my ($help, $matrixfile, $outdir, $configfile, $exptfile);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    'e|expt=s'        => \$exptfile,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrixfile, $outdir, $configfile, $exptfile){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# instantiate the object and load required data                                                          
my $cfiobj = Cfi->new;
$cfiobj->load_config ($configfile);
$cfiobj->load_aln    ($matrixfile);
my $mode = $cfiobj->get_mode;
my $root = $cfiobj->get_root;
my $charsets = $cfiobj->get_charsets;

# run expt
if ($mode eq "PARSIMONY"){
    $cfiobj->run_paup ($exptfile);
}
elsif ($mode eq "ML"){
    my @expt = split (/\//, $exptfile);
    my $exptname = pop @expt;
    my ($sample, $number, $nex) = split (/\./, $exptname);

    my $parties;
    open (F, "$exptfile");
    while (my $line = <F>){
	chomp $line;
	$parties = $line;
    }
    close (F);

    $cfiobj->sub_matrix ($sample,
			 $number,
			 $parties,
			 $matrixfile,
			 $charsets,
			 "$outdir/cmds");
    
    $cfiobj->aln_converter ("$outdir/cmds/$sample.$number.nex",
			    "nexus",
			    "$outdir/cmds/$sample.$number.phy",
			    "phylip");
    
    $cfiobj->run_raxml ("$outdir/cmds/$sample.$number.phy", 
			$root, 
			$outdir);
    
    $cfiobj->tree_converter ("$outdir/trees/RAxML_bestTree.$sample.$number.phy",
			     "newick",
			     "$outdir/trees/$sample.$number.tre",
			     "nexus");

    `rm $outdir/trees/RAxML_bestTree.$sample.$number.phy`;
    `rm $outdir/cmds/$sample.$number.nex`;
    `rm $outdir/cmds/$sample.$number.phy`;

}
else {
    print STDERR "Unknown run mode.\n";
    die;
}
