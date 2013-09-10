#!/usr/bin/perl -w

=head1 NAME

prune_tree.pl

=head1 SYNOPSIS

prune_tree.pl

Options:

--treefile is your treefile
@ARGV is a list of taxa you want to prune

The matrix must be in nexus format.

Requires the bioperl libs. 

=head1 DESCRIPTION

This program prunes taxa off a tree and outputs

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
use Bio::TreeIO;

my ($help, $treefile);
GetOptions(
    'h|help'          => \$help,
    't|tree=s'        => \$treefile,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# get the input tree
my $input = Bio::TreeIO->new(-file      => $treefile,
			     -format    => "nexus");
my $tree  = $input->next_tree();

# print the starting tree
my $treestart = Bio::TreeIO->new(-format => 'newick');
print "Starting tree:\n";
$treestart->write_tree($tree);

# prune the tree
foreach my $tax (@ARGV){
    $tree->remove_Node($tax);
}

# print the pruned tree
my $treeout = Bio::TreeIO->new(-format => 'newick');
print "Tree w/ deleted nodes:\n";
$treeout->write_tree($tree);

