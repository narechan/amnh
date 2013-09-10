#!/usr/bin/perl -w

=head1 NAME

compare_trees.pl

=head1 SYNOPSIS

compare_trees.pl

Options:

--treefile is your comparator tree
--treedir is your directory of matrix subdirs (insid)
    that you want to conmpare
--root is your root
--outdir is your output dir

Requires the bioperl libs. 

=head1 DESCRIPTION

Outputs the node frequency of all nodes found in the comparator tree
among trees in the treedir

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

my ($help, $treefile, $treedir, $root, $outdir);
GetOptions(
    'h|help'          => \$help,
    't|treefile=s'   => \$treefile,
    'd|treedir=s'   => \$treedir,
    'r|root=s'        => \$root,
    'o|outdir=s'   => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $treedir, $root, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# instantiate the object and load stuff we need
my $cfiobj = Radical->new;

# get topologies of the comparator tree
#my $topoc = $cfiobj->parse_tree ($treefile, $root, 'newick');

# load the nodes into mem                                                 
my $nodes = {};
my $treein = new Bio::TreeIO(-file => $treefile,
                             -format => 'newick');
my $treeout = new Bio::TreeIO(-file => ">$outdir/tree.index",
                              -format => 'newick');
my $tree = $treein->next_tree();

# re-root                                                                
my @ns = $tree->find_node(-id => $root);
my $n = $ns[0];
$tree->reroot($n);

# cycle through nodes and assign values                                            
open (N, ">$outdir/node.index");                                             
my $nodeindex = 0;
my @allnodes      = $tree->get_nodes;
foreach my $node (@allnodes){
    $nodeindex++;
    my $names = {};
    if ($node->is_Leaf){
        my $nodeid = $node->id;
        $names->{$nodeid} = 1;
#        if ($leaves == 1){
#            $nodes->{$nodeindex} = $names;
#            print N "$nodeindex\t$nodeid\n";
#        }
    }
    else{
        my @descendents = $node->get_all_Descendents;
        my @children    = grep{$_->is_Leaf} @descendents;

        # tips                                                             
        foreach my $leaf (@children){
            my $id = $leaf->id;
            $names->{$id} = 1;
        }
        $nodes->{$nodeindex} = $names;
        my $nodestring = join ",", sort (keys %$names);
        print N "$nodeindex\t$nodestring\n";
        $node->id($nodeindex);
    }
}
$treeout->write_tree($tree);
close (N);

# get the topologies of all the query trees
opendir (D, "$treedir");
my @subdirs = sort (readdir (D));
shift @subdirs;
shift @subdirs;
closedir (D);

my $data = {};
foreach my $subdir (@subdirs){
    print STDERR "Working on $subdir\n";

    my ($qrybases, $refbases) = split (/\-/, $subdir);

    opendir (T, "$treedir/$subdir");
    my @trees = sort (readdir (T));
    shift @trees;
    shift @trees;
    closedir (D);

    my $nodelib = {};
    foreach my $tree (@trees){
	my $topot = $cfiobj->parse_tree ("$treedir/$subdir/$tree", $root, 'nexus');
	foreach my $node (keys %$topot){
	    $nodelib->{$node}++;
	}
    }

#    open (O, ">$outdir/node.count");
    foreach my $nodec (keys %$nodes){
	my @sortedarray = sort (keys %{$nodes->{$nodec}});
	my $arrayjoin = join(",", @sortedarray);
	if ($nodelib->{$arrayjoin}){
	    $data->{$nodec}->{$qrybases}->{$refbases} = $nodelib->{$arrayjoin};
#	    print O "$nodec\t$arrayjoin\t$nodelib->{$arrayjoin}\n";
	}
	else {
	    $data->{$nodec}->{$qrybases}->{$refbases} = 0;
#	    print O "$nodec\t0\n";
	}
    }
#    close (O);
}

# print our grids for each node
foreach my $nodeindex (keys %$data){
    open (I, ">$outdir/$nodeindex.grid");
    foreach my $query (sort {$a <=> $b} keys %{$data->{$nodeindex}}){
	my @string;
	foreach my $ref (sort {$a <=> $b} keys %{$data->{$nodeindex}->{$query}}){
	    my $count = $data->{$nodeindex}->{$query}->{$ref};
	    push (@string, $count);
	}
	my $string = join "\t", @string;
	print I "$query\t$string\n";
    }
    close (I);
}
