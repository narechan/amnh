#!/usr/bin/perl -w

=head1 NAME

lild_setup.pl

=head1 SYNOPSIS

lild_setup.pl

Options:

--config       is the config file
--matrix       is the data matrix
--tree         is the total evidence tree    
--outdir       is your output dir for the run
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

my ($help, $configfile, $treefile, $matrixfile, $outdir, $window, $motion, $start, $end);
GetOptions(
    'h|help'           => \$help,
    'c|config=s'       => \$configfile,   
    't|tree=s'         => \$treefile,
    'm|matrix=s'       => \$matrixfile,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $treefile, $matrixfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}


# create dir structure for the results
`mkdir -p $outdir/cmds`;
`mkdir -p $outdir/expts`;
`mkdir -p $outdir/logs`;
`mkdir -p $outdir/res`;

#####MAIN#####

# instantiate the object and load required data
my $lildobj = Lild->new;
$lildobj->load_config ($configfile);
$lildobj->load_aln    ($matrixfile);
$lildobj->load_tree   ($treefile);

# get the nodes                                                                            
my $nodes = $lildobj->get_leaves;

# print out a nodes file                                                                            
open (LIN, ">$outdir/nodes");
foreach my $node (keys %$nodes){
    my $lineage= join (",", @{ $nodes->{$node} });
    print LIN "$node\t$lineage\n";
}
close (LIN);

# get the charsets depending on expt (slide rule or defined charsets)
my $charsets;
if ( ($window) and ($motion) and ($start) and ($end) ){
    $charsets = $lildobj->generate_partitions ($window, $motion, $start, $end);
}
else {
    $charsets = $lildobj->get_charsets;
}

# print out the expts: charset/node combos
foreach my $charset (keys %$charsets){
    foreach my $node (keys %$nodes){
	
	open (EXPT, ">$outdir/expts/$charset.$node.txt");
	my $constraintstring = join (",", @{ $nodes->{$node} });
	print EXPT "$charset\t$charsets->{$charset}\t";
	print EXPT "$node\t$constraintstring\n";
	close (EXPT);

    }
}
