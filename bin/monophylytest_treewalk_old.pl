#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::TreeIO;

my ($help, $treefile, $listfile, $root);
GetOptions(
    'h|help'          => \$help,
    't|treefile=s'   => \$treefile,
    'l|listfile=s'       => \$listfile,
    'r|root=s' => \$root,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $listfile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# store the list of taxa we're interested in
my $list = {};
open (L, "$listfile");
while (my $line = <L>){
    chomp $line;
    $list->{$line} = 1;
}
close (L);

# get the tree object and root if necessary
my $in_tree  = Bio::TreeIO->new(-file   => $treefile,
				-format => "newick");
my $tree     = $in_tree->next_tree;
if ($root){
    my @nodes = $tree->find_node(-id => $root);
    my $node = $nodes[0];
    $tree->reroot($node);
}

# get all the leaf nodes and screen for the ones
# with our taxa of interest
my @seqskept;
my @leaveskept;
my @leaves = $tree->get_leaf_nodes;
foreach my $leaf (@leaves){
    my $id = $leaf->id;
    my ($taxon, $acc) = split (/\#/, $id);
    if (exists ($list->{$taxon})){
	push (@seqskept, $id);
	push (@leaveskept, $leaf);
    }
    else{
	next;
    }
}

# for every leaf walk up the tree until we hit a node
# where the opposite lineage contains something other than 
# our taxa of interest
my $clades = {};
foreach my $leaf (@leaveskept){
    my $leafid = $leaf->id;

    # get the leaf's path all the down to the root
    my @path;
    while (defined $leaf){
	push (@path, $leaf);
	$leaf = $leaf->ancestor;
    }

    # take the leaf off the path
    my $pathcache = shift @path;
    
    # query the descendants of each node up the path
    my $largestnode;
    foreach my $node (@path){
	my @descendents = $node->get_all_Descendents;
	my @children = grep{$_->is_Leaf} @descendents;
	
	# check to see if the ids for this node
	# contain only the desired taxa
	my $signal = 0;
	foreach my $leaf (@children){
	    my $id = $leaf->id;
	    my ($taxon, $acc) = split (/\#/, $id);
	    
	    # check to see if the id is the list of desired taxa
	    if (exists ($list->{$taxon})){
		next;
	    }
	    else {
		$signal = 1;
		last;
	    }
	}
	
	if ($signal == 1){
	    $largestnode = $pathcache;
	    last;
	}
	else {
	    $pathcache = $node;
	}
    }
    
    # get all the taxa in the largest node and record
    my @largedescendents = $largestnode->get_all_Descendents;
    my @largechildren = grep{$_->is_Leaf} @largedescendents;
    my @llseqs;

    # if there are no descendants, this is a singleton clade
    # otherwise print all leaves to define the clade
    my $numchildren = @largechildren;
    if ($numchildren == 0){
	$clades->{$leafid} = 1;
#	print "$leafid\t$leafid\n";
    }
    else{
	foreach my $largeleaf (@largechildren){
	    my $llid = $largeleaf->id;
	    push (@llseqs, $llid);
	}
	my $llseqs = join ",", sort @llseqs;
	$clades->{$llseqs} = 1;
#	print "$leafid\t$llseqs\n";
    }
}
	
# print all the unique largest clades that 
# obey the largest monophyly rules
my $cladecounter = 0;
foreach my $clade (keys %$clades){
    $cladecounter++;
    print "clade$cladecounter\t$clade\n";
}
