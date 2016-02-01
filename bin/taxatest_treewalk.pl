#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::TreeIO;

my ($help, $treefile, $qrylistfile, $tglistfile, $root);
GetOptions(
    'h|help'          => \$help,
    't|treefile=s'   => \$treefile,
    'a|qrylistfile=s'       => \$qrylistfile,
    'b|tglistfile=s' => \$tglistfile,
    'r|root=s' => \$root,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $qrylistfile, $tglistfile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# store the list of taxa we're interested in targeting
my $tglist = {};
open (L, "$tglistfile");
while (my $line = <L>){
    chomp $line;
    $tglist->{$line} = 1;
}
close (L);

# also store the qrygroups
my $qrylist = {};
open (I, "$qrylistfile");
while (my $line = <I>){
    chomp $line;
    $qrylist->{$line} = 1;
}
close (I);

# get the tree object and root if necessary
my $in_tree  = Bio::TreeIO->new(-file   => $treefile,
				-format => "newick");
my $tree     = $in_tree->next_tree;
if ($root){
    my @nodes = $tree->find_node(-id => $root);
    my $node = $nodes[0];
    $tree->reroot($node);
}

# get all the qrygrp leaf nodes
my @seqskept;
my @leaveskept;
my @leaves = $tree->get_leaf_nodes;
foreach my $leaf (@leaves){
    my $id = $leaf->id;
    my ($taxon, $acc) = split (/\#/, $id);
    if (exists ($qrylist->{$taxon})){
	push (@seqskept, $id);
	push (@leaveskept, $leaf);
    }
    else{
	next;
    }
}

# for every qry leaf walk up the tree until we hit a node
# with some target taxa, and record
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
    my $signal = 0;
    my @targetseqs;
    foreach my $node (@path){
	my @descendents = $node->get_all_Descendents;
	my @children = grep{$_->is_Leaf} @descendents;
	
	# check to see if the ids for this node
	# contain any target sequences
	foreach my $leaf (@children){
	    my $id = $leaf->id;
	    my ($taxon, $acc) = split (/\#/, $id);
	    
	    # if you haven't yet encountered an ingroup,
	    # keep searching for outgroups
	    if (exists ($tglist->{$taxon})){
		$signal++;
		push (@targetseqs, $id);
	    }
	    else{
		next;
	    }
	
	}
	
	if ($signal > 0){
	    last;
	}
	else {
	    next;
	}
    }
    
    my $targetseqs = join ";", @targetseqs;
    print "$treefile\t$leafid\t$targetseqs\n";
}
