#!/usr/bin/perl -w

# -i is the input tree file
# -f is the input format
# -r is the root (optional)

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::TreeIO;

my ($infile, $informat, $root);
GetOptions(
	   'i|infile=s'      => \$infile,
	   'f|informat=s'    => \$informat,
	   'r|root=s'        => \$root,
	   );

=head
my $treeio = Bio::TreeIO->new(-format => "$informat",
			      -file => "$infile");

while( my $tree = $treeio->next_tree ) {
    for my $node ( $tree->get_nodes ) {
	my $nodeid = $node->id;
	my $boot   = $node->bootstrap;
	print "$nodeid\t$boot\n";
    }
}
=cut

my $treein = new Bio::TreeIO(-file => "$infile",
			     -format => "$informat");
my $tree = $treein->next_tree();

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
my $internal_counter = 0;
foreach my $node (@internalnodes){
    $internal_counter++;

    my $nodeid = $node->id; #these are the bootstraps for all internal nodes
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
    my $sortedarrayCnt = @sortedarray;
    my $arrayjoin = join(",", @sortedarray);
    $lineages->{$arrayjoin} = [@sortedarray];
    
    if ($nodeid){
	print "$internal_counter\t$arrayjoin\t$nodeid\t$sortedarrayCnt\n";
    }
    else {
	print "$internal_counter\t$arrayjoin\tNA\t$sortedarrayCnt\n";
    }
}
