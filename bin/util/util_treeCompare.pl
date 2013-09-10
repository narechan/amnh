#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Cfi;

my ($help, $TEtree, $treedir, $root);
GetOptions(
    'h|help'          => \$help,
    't|TEtree=s'      => \$TEtree,
    'd|treedir=s'     => \$treedir,
    'r|root=s'        => \$root,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($TEtree, $treedir, $root){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# instantiate the object and load required data                                                           
my $cfiobj = Cfi->new;
my $topoTE  = $cfiobj->parse_tree ($TEtree, $root);

opendir (TD, "$treedir");
my @trees = sort (readdir (TD));
shift @trees;
shift @trees;
closedir (TD);

# generate cfi data, and do topo comparisons
my $TEstats = {};
foreach my $tree (@trees){
    print STDERR "$tree\n";

    $cfiobj->tree_converter ("$treedir/$tree", "newick", "/tmp/treefile.apurva", "nexus");
    my $cfitopo = $cfiobj->parse_tree ("/tmp/treefile.apurva", $root);
    my ($t, $topocomp) = $cfiobj->compare_all_nodes ($topoTE, $cfitopo);

    foreach my $top (sort keys %$t){
	($TEstats->{$top}++) if ($t->{$top} == 1);
    }
    `rm /tmp/treefile.apurva`;
}

foreach my $top (sort keys %$TEstats){
    print "$top\t$TEstats->{$top}\n";
}
