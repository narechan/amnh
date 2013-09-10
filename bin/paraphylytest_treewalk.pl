#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::TreeIO;

my ($help, $treefile, $outlistfile, $inlistfile, $root);
GetOptions(
    'h|help'          => \$help,
    't|treefile=s'   => \$treefile,
    'a|outlistfile=s'       => \$outlistfile,
    'b|inlistfile=s' => \$inlistfile,
    'r|root=s' => \$root,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $outlistfile, $inlistfile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# store the list of taxa we're interested in
# (in this case the list is a list of outgroups)
my $outlist = {};
open (L, "$outlistfile");
while (my $line = <L>){
    chomp $line;
    $outlist->{$line} = 1;
}
close (L);

# also store the ingroups
my $inlist = {};
open (I, "$inlistfile");
while (my $line = <I>){
    chomp $line;
    $inlist->{$line} = 1;
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

# get all the outgrp leaf nodes and screen for the ones
# with our taxa of interest
my @seqskept;
my @leaveskept;
my @leaves = $tree->get_leaf_nodes;
foreach my $leaf (@leaves){
    my $id = $leaf->id;
    my ($taxon, $acc) = split (/\#/, $id);
    if (exists ($outlist->{$taxon})){
	push (@seqskept, $id);
	push (@leaveskept, $leaf);
    }
    else{
	next;
    }
}

# for every leaf walk up the tree until we hit a node
# with some ingroup taxa, then continue walking up until
# another outgroup is encountered. stop there and record 
# as paraphyletic group
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
	# contain only outgroups
	my $signal = 0;
	foreach my $leaf (@children){
	    my $id = $leaf->id;
	    my ($taxon, $acc) = split (/\#/, $id);
	    
	    # if you haven't yet encountered an ingroup,
	    # keep searching for outgroups
	    if ($signal == 0){
		if (exists ($outlist->{$taxon})){
		    next;
		}
		else {
		    $signal++;
		    last;
		}
	    }

	    # if you have encountered an ingroup,
	    # keep searching for ingroups
	    elsif ($signal == 1){
		if (exists ($inlist->{$taxon})){
                    next;
                }
                else {
                    $signal++;
		    last;
                }
	    }
	    else {
		next;
	    }
	}
	
	if ($signal == 2){
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
