#!/usr/bin/perl -w

=head1 NAME

treeSupports_setup.pl

=head1 SYNOPSIS

treeSupports_setup.pl

Options:

--config       is the config file
--matrix       is the data matrix
--tree         is the total evidence tree    
--outdir       is your output dir for the run
--window       is your window size if doing                                       
    sildeRule (optional)                                         
--motion       is the number of residues the                                 
    window moves per calculation (5' -> 3')                           
    if doing slideRule (optional)
--start        is the start position of slideRule
    (optional)
--end          is the end position of slideRule
    (optional)

Both the matrix and the tree must be in
nexus format.

Dependencies:                                                                                        
                                                                                               
PAUP must be installed and in $PATH    
Requires the bioperl libs.
Requires TreeSupports.pm

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2012 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use TreeSupports;

my ($help, $treefile, $matrixfile, $configfile, $outdir, $window, $motion, $start, $end);
GetOptions(
    'h|help'           => \$help,
    't|tree=s'         => \$treefile,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $matrixfile, $configfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# instantiate the object and load required data
my $supportobj = TreeSupports->new;
$supportobj->load_config ($configfile);
my $root = $supportobj->get_root;
my $config = $supportobj->get_config;

my %tests;
my @tests = split (/\s/, $config->{'TESTS'});
foreach my $test (@tests){
    $tests{$test} = 1;
}

$supportobj->load_aln    ($matrixfile);
$supportobj->load_tree   ($treefile, $root);
my $nodes    = $supportobj->get_leaves;
my $lineages = $supportobj->get_lineages;

# print out a nodes file
open (N, ">$outdir/nodes");
foreach my $node (sort keys %$nodes){
    my $lineage = join (",", @{ $nodes->{$node} });
    print N "$node\t$lineage\n";
}
close (N);

# print out a lineages file
open (L, ">$outdir/lineages");
foreach my $node (sort keys %$lineages){
    my @lins;
    foreach my $lin (sort keys %{ $lineages->{$node} }){
        push (@lins, join (",", @{ $lineages->{$node}->{$lin} }));
    }
    my $nodestream = join ("\t", @lins);
    print L "$node\t$nodestream\n";
}
close (L);

# get the charsets depending on the expt (slide rule or defined char sets)
# if slide rule, those charsets override those defined in the nxs file
my $charsets;
if ( ($window) and ($motion) and ($start) and ($end) ){
    $charsets = $supportobj->generate_partitions ($window, $motion, $start, $end);
}
else {
    $charsets = $supportobj->get_charsets;
}

# NDI expts
if (exists ($tests{'ndi'})){
    `mkdir -p $outdir/ndi/expts`;
    `mkdir -p $outdir/ndi/cmds`;
    `mkdir -p $outdir/ndi/logs`;
    `mkdir -p $outdir/ndi/trees`;
    foreach my $node (keys %$nodes){
	open (NDI, ">$outdir/ndi/expts/$node.txt");
	my $constraintstring = join (",", @{ $nodes->{$node} });
	print NDI "$node\t$constraintstring\n";
	close (NDI);
    }
}

# BS expts
if (exists ($tests{'bs'})){
    `mkdir -p $outdir/bs/expts`;
    `mkdir -p $outdir/bs/cmds`;
    `mkdir -p $outdir/bs/logs`;
    `mkdir -p $outdir/bs/prunes`;
    `mkdir -p $outdir/bs/trees`;
    foreach my $charset (keys %$charsets){
	open (BS, ">$outdir/bs/expts/$charset.txt");
	print BS "$charset\t$charsets->{$charset}\n";
	close (BS);
    }
}
    
# PBS expts
if (exists ($tests{'pbs'})){
    `mkdir -p $outdir/pbs/expts`;
    `mkdir -p $outdir/pbs/cmds`;
    `mkdir -p $outdir/pbs/logs`;
    `mkdir -p $outdir/pbs/trees`;
    foreach my $node (keys %$nodes){
	open (PBS, ">$outdir/pbs/expts/$node.txt");
	my $constraintstring = join (",", @{ $nodes->{$node} });
	print PBS "$node\t$constraintstring\n";
	close (PBS);
    }
}

# unconstrained PBS expt
open (UNCON, ">$outdir/pbs/expts/unconstrained.txt");
print UNCON "unconstrained\n";
close (UNCON);
