#!/usr/bin/perl -w

=head1 NAME

compare_trees_bio_phylo.pl

=head1 SYNOPSIS

compare_trees_bio_phylo.pl

Options:

--treefile1 is your first treefile
--treefile2 is your second treefile

--treedir is your directory of trees

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
use Bio::Phylo::IO 'parse';
use Data::Dumper;

my ($help, $treefile1, $treefile2, $treedir);
GetOptions(
    'h|help'          => \$help,
    'a|treefile1=s'   => \$treefile1,
    'b|treefile2=s'   => \$treefile2,
    't|treedir=s'     => \$treedir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile1, $treefile2, $treedir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# get trees
opendir (TD, "$treedir");
my @trees = sort (readdir (TD));
shift @trees;
shift @trees;
closedir (TD);

my $trees = {};
my $ladtrees = {};
foreach my $treefile (@trees){
#    print STDERR "$treefile\n";
    my $proj = parse(
		     -file => "$treedir/$treefile",
		     -format => "nexus",
		     -as_project => 1,
		     );
    
    my $forest = $proj->get_forests->[0];
    my $tree   = $forest->get_entities->[0];
    
#    my $root = $tree->get_root->get_name;

#    my $treest  = $tree->to_newick;
#    $trees->{$treest}++;
#    push @{$trees->{$treest}}, $treefile;

    my $ladtree = $tree->ladderize;
    my $ladtreest = $ladtree->to_newick;
#    $ladtrees->{$ladtreest}++;
#    push @{$ladtrees->{$ladtreest}}, $treefile;

    print "$treefile $ladtreest\n";
}

my $tcounter = 0;
foreach my $tree (keys %$trees){
    $tcounter++;
    my @treearray = @{$trees->{$tree}};
    my $treearray = @treearray;
    my $treearrayst = join ",", @treearray;
    print "$tcounter\t$tree\t$treearray\t$treearrayst\n";
}

print "---------\n";

my $ltcounter = 0;
foreach my $tree (keys %$ladtrees){
    $ltcounter++;
    my @treearray = @{$ladtrees->{$tree}};
    my $treearray = @treearray;
    my $treearrayst = join ",", @treearray;
    print "$ltcounter\t$tree\t$treearray\t$treearrayst\n";
}

sub parse_tree{
    my $treefile1 = shift;

    my $proj = parse(
                     -file => $treefile1,
                     -format => "nexus",
                     -as_project => 1,
                     );

    my $forest = $proj->get_forests->[0];
    my $tree   = $forest->get_entities->[0];

    my $lineages = {};
    for my $node (@{$tree->get_internals}){
        my @terminals = @{$node->get_terminals};

        my @array;
        foreach my $term (@terminals){
            push (@array, $term->get_name);
        }

        my @sortedarray = sort @array;
        my $arrayjoin = join(",", @sortedarray);
        $lineages->{$arrayjoin} = [@sortedarray];
    }
    return ($lineages);
}

