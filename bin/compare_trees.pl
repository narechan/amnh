#!/usr/bin/perl -w

=head1 NAME

compare_trees.pl

=head1 SYNOPSIS

compare_trees.pl

Options:

--treefile1 is your first treefile
--treefile2 is your second treefile
--root is your root

The matrix must be in nexus format.

Requires the bioperl libs. 

=head1 DESCRIPTION

This computes the consensus fork between two trees

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
use Radical;
use Data::Dumper;

my ($help, $treefile1, $treefile2, $root);
GetOptions(
    'h|help'          => \$help,
    'a|treefile1=s'   => \$treefile1,
    'b|treefile2=s'   => \$treefile2,
    'r|root=s'        => \$root,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile1, $treefile2, $root){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $cfiobj = Radical->new;

# get topologies of the pruned trees
my $topo1 = $cfiobj->parse_tree ($treefile1, $root, 'nexus');
my $topo2 = $cfiobj->parse_tree ($treefile2, $root, 'nexus');

# check to see if the topologies of the pruned trees are the same
my $cfidata = 0;
my ($t, $topocomp) = $cfiobj->compare_all_nodes ($topo1, $topo2);
print  Dumper ($t);
my $counter = 0;
foreach my $top (sort keys %$t){
    $counter++;
    $cfidata += $t->{$top};
}
$cfidata--; 
$cfidata--;

print "$topocomp\t$cfidata\n";
print "$counter\n";
