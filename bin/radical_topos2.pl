#!/usr/bin/perl -w

=head1 NAME

radical_topos.pl

=head1 SYNOPSIS

  radical_topos.pl -- Mine the RADICAL analysis for topologies.

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --outdir      Is your output dir
 --config      Is the configuration for your tests

The aln file must be in nexus format.

Dependencies:

Requires the bioperl libraries
Requires Radical.pm

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2011 American Museum of Natural History

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
use Bio::Phylo::IO 'parse';

my ($help, $matrixfile, $outdir, $configfile);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrixfile, $outdir, $configfile){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# instantiate the object and load required data                                                           
my $cfiobj = Radical->new;
$cfiobj->load_config ($configfile);
$cfiobj->load_aln    ($matrixfile);

# get the root from config
my $root    = $cfiobj->get_root;
my $sections = $cfiobj->get_sections;

# get charset number
my $charsets = $cfiobj->get_charsets;
my @charsets = keys (%$charsets);
my $charsetnum = @charsets;

# get all the trees                                                                                  
opendir (TREES, "$outdir/trees");
my @trees = sort (readdir (TREES));
shift @trees;
shift @trees;
closedir (TREES);

# find unique trees
my $trees = {};
foreach my $treefile (@trees){
    print STDERR "$treefile\n";
    my $proj = parse(
                     -file => "$outdir/trees/$treefile",
                     -format => "nexus",
                     -as_project => 1,
                     );

    my ($forest) = @{$proj->get_forests};
    my ($tree)   = @{$forest->get_entities};

    my $ladtree = $tree->ladderize;
    my $ladtreest = $ladtree->to_newick;
#    push @{$trees->{$ladtreest}}, $treefile;
    $trees->{$ladtreest}++;
}

# print the unique trees and their freqs
my $tcounter = 0;
foreach my $tree (keys %$trees){
    $tcounter++;
#    my @treearray = @{$trees->{$tree}};
#    my $treearray = @treearray;
#    my $treearrayst = join ",", @treearray;
#    print "$tcounter\t$tree\t$treearray\n"; #$treearrayst\n";
    print "$tcounter\t$tree\t$trees->{$tree}\n";
}

# build the consensus tree across all reps and additions                                                 
print STDERR "Consensus Building\n";
$cfiobj->generate_consensus ("$outdir", $matrixfile, \@trees, "all");
$cfiobj->run_paup ("$outdir/consensus/all.nex");

# bin the trees
my $treesects = {};
my $start = 1;
my $end = 0;
if ($sections =~m/\,/){
    $sections =~s/\s//g;
    my @sections = split (/\,/, $sections);
    foreach my $sect (@sections){
        my $sestring = $start . "-" . $sect;
        $treesects->{$sestring} = [];
        $start = $sect + 1;
    }
}
else {
    my $sectparts = round ($charsetnum / $sections);

    until ($start >= $charsetnum){
        $end = $start + $sectparts - 1;
        if ($end > $charsetnum){
            $end = $charsetnum;
        }

        my $sestring = $start . "-" . $end;
        $treesects->{$sestring} = [];
        $start = $end + 1;
    }
}

foreach my $tree (@trees){
    my ($sample, $part, $tre) = split (/\./, $tree);

    foreach my $sect (sort keys %$treesects){
        my ($s, $e) = split (/\-/, $sect);
        if (($part >= $s) and ($part <= $e)){
            push @{$treesects->{$sect}}, $tree;
        }
        else {
            next;
        }
    }
}

# build the consensus in each bin                                                                       
# and harvest all nodes and all nodes in each bin                                                          
my $allnodes = {};
my $secnodes = {};
foreach my $sect (sort keys %$treesects){
    $cfiobj->generate_consensus ("$outdir", $matrixfile, \@{$treesects->{$sect}}, $sect);
    $cfiobj->run_paup ("$outdir/consensus/$sect.nex");

    foreach my $tree (@{$treesects->{$sect}}){
	my $cfitopo = $cfiobj->parse_tree ("$outdir/trees/$tree", $root);
        foreach my $n (sort keys %$cfitopo){
            $allnodes->{$n}++;
            $secnodes->{$sect}->{$n}++;
        }
    }
}

# print all harvested nodes                                                                               
print STDERR "Node distributon\n";
open (AN, ">$outdir/summaries/all.nodes");
foreach my $node (sort {$allnodes->{$b} <=> $allnodes->{$a}} keys %$allnodes){
    print AN "$node\t$allnodes->{$node}\n";
}
close (AN);

# print all harvested nodes in the curve sections                                                       
foreach my $sec (sort keys %$secnodes){
    open (SN, ">$outdir/summaries/$sec.nodes");
    foreach my $node (sort {$secnodes->{$sec}->{$b} <=> $secnodes->{$sec}->{$a}} keys %{$secnodes->{$sec}}){
        print SN "$node\t$secnodes->{$sec}->{$node}\n";
    }
    close (SN);
}

####SUBS####

sub round {
    my($number) = shift;
    return int($number + .5);
}
