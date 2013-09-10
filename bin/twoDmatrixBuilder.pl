#!/usr/bin/perl -w

=head1 NAME

twoDmatrixBuilder.pl

=head1 SYNOPSIS
    
  twoDmatrixBuilder.pl -- the program will partition a presence absence matrix according to 
    taxonomic node from the phylogeny. Presence/absence data in this case is location data.

Options:

 --help        Show brief help and exit
 --treefile    Taxonomic treefile   
 --matrix      Presence absence matrix
 --lookup      File the coordinates PA position with taxon
 --outdir      Output dir

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

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::TreeIO;

my ($help, $treefile, $matrix, $lookup, $outdir);
GetOptions(
    'h|help'          => \$help,
    't|treefile=s'    => \$treefile,
    'm|matrix=s'      => \$matrix,
    'l|lookup=s'      => \$lookup,
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;

pod2usage if $help;

for my $option ($treefile, $matrix, $lookup, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# store data in the lookup file
# matrix position to taxon
my $table = {};
open (LU, "$lookup");
while (my $line = <LU>){
    chomp $line;
    
    my ($pos, $acc) = split (/\t/, $line);
    $table->{$acc} = $pos;
}
close (LU);

# get the tree topology and build nodal partitions
my $charsets = {};
my $counter = 0;
my ($topo, $ultrametrics) = parse_tree ("$treefile");

open (NDS, ">$outdir/nodes");
open (LNS, ">$outdir/lengths");
foreach my $node (sort keys %$topo){
    $counter++;
    print NDS "$counter\t$node\n";
    print LNS "$counter\t$ultrametrics->{$node}\n";

    foreach my $leaf (@{$topo->{$node}}){
	push @{$charsets->{$counter}}, $table->{$leaf};
    }
}
close (NDS);
close (LNS);

# build the new nexus file	
`cp $matrix $outdir`;
open (PRTS, ">$outdir/partitions");
print PRTS "Begin Sets;\n";
foreach my $set (sort {$a <=> $b} %$charsets){
    next unless ($set =~m/^[+-]?\d+$/); #HACK -- for weird value inserted into $charsets

    print STDERR "$set\n";
    my $charstring = join " ", @{$charsets->{$set}};
    print PRTS "Charset Part$set = $charstring;\n";
}
print PRTS "End;\n";
close (PRTS);

`cat $outdir/$matrix $outdir/partitions > $outdir/nexus`;

###SUBS###

sub parse_tree {
    my $treefile = shift;
    my $root     = shift;

    my $input = Bio::TreeIO->new(-file      => $treefile,
                                 -format    => "nexus");

    # now remove the taxa and analyze the culled trees                                            
    my $tree = $input->next_tree();

    # re-root for nodal comparisons, if required                                                    
    if ($root){
        my @nodes = $tree->find_node(-id => $root);
        my $node = $nodes[0];
	$tree->reroot($node);
    }

    # get all nodes and internal nodes                                                             
    my @allnodes      = $tree->get_nodes;
    my @internalnodes = grep {!$_->is_Leaf} @allnodes;

    my $lineages = {};
    my $blengths = {};
    my $internal_counter = 0;
    foreach my $node (@internalnodes){
        $internal_counter++;
	
        my @descendents = $node->get_all_Descendents;
	my @children    = grep{$_->is_Leaf} @descendents;

        my @array;
        foreach my $leaf (@children){
            my $id = $leaf->id;
            push (@array, $id);
        }

        # this checks to see that all nodes have a uniq                                                
	# list of children --> compensates for a pruning bug                                          
        # in TreeIO where an extra set of parens is left                                             
	# on outermost nest;                        
	my @sortedarray = sort @array;
        my $arrayjoin = join(",", @sortedarray);
        $lineages->{$arrayjoin} = [@sortedarray];

#	my $blength = $node->branch_length;
	my $blength = $node->height;
	$blengths->{$arrayjoin} = $blength;
    }

    return ($lineages, $blengths);
}
