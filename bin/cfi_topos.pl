#!/usr/bin/perl -w

=head1 NAME

cfi_topos.pl

=head1 SYNOPSIS

  cfi_topos.pl

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --outdir      Is your output dir
 --config      Is the configuration for your tests

The aln file must be in nexus format.

PAUP must be in your path
Requires Cfi.pm

The config must specify all parameters for the
search (treecommand), the reps (samples), 
the max trees searched (maxtrees).

Maxtrees can be left undefined in which case it will be
set to 100 (default) and the command, "SET increase=auto" 
will be used.
Examp:

TREECOMMAND=hsearch swap=nni addseq=closest
MAXTREES=200
SAMPLES=10

=head1 DESCRIPTION

Calculate CFI given successive additions of all partitions,
one at a time

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2009 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Cfi;
use Statistics::Descriptive;

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
my $cfiobj = Cfi->new;
$cfiobj->load_config ($configfile);
$cfiobj->load_aln    ($matrixfile);

# get the root from config
my $root    = $cfiobj->get_root;
my $sections = $cfiobj->get_sections;
#my $mode = $cfiobj->get_mode;

#my $treetype;
#if ($mode eq "PARSIMONY"){
#    $treetype = "nexus";
#}
#elsif ($mode eq "ML"){
#    $treetype = "newick";
#}
#else {
#    print STDERR "Unknown analysis type.\n";
#    die;
#}

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

# create a distribution of all tree topologies                                                         
print STDERR "Topology Distribution\n";
my $topodistfiles = {};
my $topodisttopos = {};
my $topodistcounter = 1;
my $topocatfiles = {};
my $topocattopos = {};
my $topocatcounter = 1;

# evaulate the first tree to seed the comparisons
my ($rep, $pos, $suf) = split (/\./, $trees[0]);
my $firsttopo = $cfiobj->parse_tree ("$outdir/trees/$trees[0]", $root);

$topodisttopos->{$topodistcounter} = $firsttopo;
push @{$topodistfiles->{$topodistcounter}}, $trees[0];

$topocattopos->{$pos}->{$topocatcounter} = $firsttopo;
push @{$topocatfiles->{$pos}->{$topocatcounter}}, $trees[0];

# cycle through all remaining trees
foreach my $tree (@trees){
    print STDERR "$tree\n";
    next if ($tree eq $trees[0]);

    my ($rep, $pos, $suf) = split (/\./, $tree);

    my $topo = $cfiobj->parse_tree ("$outdir/trees/$tree", $root);
    my $signal = 0;

    # cycle through all unique topologies saved
    foreach my $savedtopo (keys %$topodisttopos){
        my ($t, $topocomp) = $cfiobj->compare_all_nodes ($topo, $topodisttopos->{$savedtopo});

	# if the trees are the same, do the following
        if ($topocomp == 1){
            push @{$topodistfiles->{$savedtopo}}, $tree;
            $signal = 1;
            last;
        }
        else {
            next;
        }
    }

    # if the trees are different, do the following
    if ($signal == 0){
        $topodistcounter++;
        $topodisttopos->{$topodistcounter} = $topo;
        push @{$topodistfiles->{$topodistcounter}}, $tree;
    }

    my $catsignal = 0;
    
    # cycle through the unique topologies saved for this pos
    foreach my $savedcattopo (keys %{$topocattopos->{$pos}}){
	my ($t, $topocomp) = $cfiobj->compare_all_nodes ($topo, $topocattopos->{$pos}->{$savedcattopo});
	
        # if the trees are the same, do the following                                                     
        if ($topocomp == 1){
            push @{$topocatfiles->{$pos}->{$savedcattopo}}, $tree;
            $catsignal = 1;
            last;
        }
        else {
            next;
        }
    }

    # if the trees are different, do the following                                                    
    if ($catsignal == 0){
        $topocatcounter++;
        $topocattopos->{$pos}->{$topocatcounter} = $topo;
        push @{$topocatfiles->{$pos}->{$topocatcounter}}, $tree;
    }
    
}

# print out the overall tree topos
open (TD, ">$outdir/topology.distribution");
foreach my $topoindex (sort keys %$topodisttopos){
    print TD "$topoindex\t";

    my @branches = sort (keys %{$topodisttopos->{$topoindex}});
    my $topology = join ";", @branches;
    print TD "$topology\t";

    my $filecount = @{$topodistfiles->{$topoindex}};
    print TD "$filecount\n";
}
close (TD);

# print out the position specific tree topo distribution
open (PD, ">$outdir/topology_concat.distribution");
foreach my $catindex (sort {$a <=> $b} keys %$topocatfiles){
    my @catuniqs = keys %{$topocatfiles->{$catindex}};
    my $catuniqs = @catuniqs;
    print PD "$catindex\t$catuniqs\n";
}
close (PD);

=head
# build the consensus tree across all reps and additions                                            
print STDERR "Consensus Building\n";
$cfiobj->generate_consensus ("$outdir", $matrixfile, \@trees, "all");
$cfiobj->run_paup ("$outdir/consensus/all.nex");

# build consensus trees for each curve partition                                                    
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

# bin the trees                                                                                            
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
=cut
sub round {
    my($number) = shift;
    return int($number + .5);
}
