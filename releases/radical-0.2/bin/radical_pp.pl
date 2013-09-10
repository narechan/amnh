#!/usr/bin/perl -w

=head1 NAME

radical_pp.pl

=head1 SYNOPSIS

  radical_pp.pl -- Post process the RADICAL tree libraries to generate stats
    consensus trees, and radical curves.

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --outdir      Is your output dir
 --config      Is the configuration for your tests
 --purge       If you want to purge the working files to save disk space (optional)
 --tree        Is your TE treefile (optional)                                                               
 --dist        Is your topology distribution (optional)
 --index       Is the index to analyze in the topo distribution (optional)     

Note that either --tree must be set or --index/--dist must be set.
The aln file must be in nexus format.

Dependencies:

Requires the bioperl libraries
Requires Radical.pm
Requires Statistics::Descriptive
Requires Bio::Phylo

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
use Statistics::Descriptive;
use Bio::Phylo::IO 'parse';

my ($help, $matrixfile, $treefile, $outdir, $configfile, $purge, $index, $dist);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    't|tree=s'        => \$treefile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    'p|purge'         => \$purge,
    'i|index=s'         => \$index,	   
    'd|dist=s'          => \$dist,
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

# get charset number and ntax
my $ntax     = $cfiobj->get_ntax;
my $nodes    = $ntax - 3;
my $charsets = $cfiobj->get_charsets;
my @charsets = keys (%$charsets);
my $charsetnum = @charsets;
my $sections = $cfiobj->get_sections;

# get TE topology and the root from config
my $root    = $cfiobj->get_root;

($index = "$treefile") unless ($index);
my $topoindex = {};
if ($treefile){
    $topoindex  = $cfiobj->parse_tree ($treefile, $root, "nexus");
}
elsif ($dist and $index){
    open (TD, "$dist");
    while (my $line = <TD>){
	chomp $line;
	my ($num, $topo) = split (/\t/, $line);
	
	if ($num == $index){
	    open (NTF, ">$outdir/$index.newick.tre");
	    print NTF "$topo;\n";
	    close (NTF);
	    $topoindex  = $cfiobj->parse_tree ("$outdir/$index.newick.tre", $root, "newick");
	    last;
	}
	else {
	    next;
	}
    }
    close (TD);
}
else {
    print STDERR "Need a topology!\n";
    die;
}

# make holding directories
`mkdir -p $outdir/summaries/$index`;
`mkdir -p $outdir/output/$index/samples`;

# get trees                                                                                
opendir (L, "$outdir/trees");
my @trees = sort (readdir(L));
shift @trees;
shift @trees;
closedir (L);

# build consensus trees and tabulate nodes unless already done 
# for this radical library (avoids redundancy and file conflict in 
# forked processes)
unless (-e "$outdir/consensus/all.nex"){
    
    # build the consensus tree across all reps and additions                
    print STDERR "Consensus Building and node distribution\n";
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
	next unless (@{$treesects->{$sect}});
	$cfiobj->generate_consensus ("$outdir", $matrixfile, \@{$treesects->{$sect}}, $sect);
	$cfiobj->run_paup ("$outdir/consensus/$sect.nex");
	
	foreach my $tree (@{$treesects->{$sect}}){
	    my $cfitopo = $cfiobj->parse_tree ("$outdir/trees/$tree", $root, "nexus");
	    foreach my $n (sort keys %$cfitopo){
		$allnodes->{$n}++;
		$secnodes->{$sect}->{$n}++;
	    }
	}
    }
    
    # print all harvested nodes                                                    
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
}

# generate cfi data, and do topo comparisons
print STDERR "$index: Doing CFI and storing data\n";
my $cfi = {};
my $topo = {};
my $cfitv = {};
my $topotv = {};
my $partnum = 0;
foreach my $tree (@trees){
    my ($sample, $included, $nex) = split (/\./, $tree, 3);

    # calculate the number of concats done (max concat pt)
    ($partnum = $included) if ($included > $partnum);
    
    my $cfitopo = $cfiobj->parse_tree ("$outdir/trees/$tree", $root, "nexus");
    my ($t, $topocomp) = $cfiobj->compare_all_nodes ($topoindex, $cfitopo);
    
    my $cfidata = 0;
    foreach my $top (sort keys %$t){
	push @{$topo->{$top}->{$included}}, $t->{$top};
	$topotv->{$top}->{$sample}->{$included} = $t->{$top};
	$cfidata += $t->{$top};
    }
    $cfidata--; 
    $cfidata--;
    push @{$cfi->{$included}}, $cfidata;
    $cfitv->{$sample}->{$included} = $cfidata;
}
$partnum--;

# create summaries for straight cfi
# and get cfi fixation point
print STDERR "Sum CFI\n";
my $fixation = {};
my $asymcfi = 
    summarize ($cfi, $cfitv, $partnum, $nodes, 0, $outdir, $index); #zero is the cfi node
$fixation->{'CFI'} = $asymcfi;

# create summaries for each node in the topology
# and create a node reference file
# and get node fixation points
my $counter = 0;
my $topocode = {};
open (NODES, ">$outdir/summaries/$index/nodes");
foreach my $node (sort keys %$topo){
    print STDERR "$index: Sum $node\n";
    $counter++;
    $topocode->{$node} = $counter;
    print NODES "$counter\t$node\n";

    my $asymnode = 
	summarize ($topo->{$node}, $topotv->{$node}, $partnum, 1, $counter, $outdir, $index);
    $fixation->{$node} = $asymnode;
}
close (NODES);

# print out node fixation points table
open (FIX, ">$outdir/summaries/$index/fixation.points");
foreach my $nd (sort keys %$fixation){
    my $asymstring = join "\t", @{$fixation->{$nd}};
    print FIX "$nd\t$asymstring\n";
}
close (FIX);

# create the final summary files
opendir (D, "$outdir/output/$index");
my @files = sort (readdir (D));
shift @files;
shift @files;
closedir (D);

# transform
my $datafin = {};
foreach my $file (@files){
    next if ($file eq "samples");
    
    my ($k, $n) = split (/\./, $file);
    
    open (F, "$outdir/output/$index/$file");
    if ($k eq "auc"){
	while (my $line = <F>){
	    chomp $line;
	    my ($samp, $wt, $loess) = split (/\t/, $line);
	    push @{$datafin->{$k}->{$samp}->{$n}}, $wt, $loess;
	}
    }
    elsif (($k eq "coords") or ($k eq "coords_loess")){
	while (my $line = <F>){
	    chomp $line;
	    my ($part, $cfi) = split (/\t/, $line);
	    $datafin->{$k}->{$part}->{$n} = $cfi;
	}
    }
    elsif ($k eq "integration"){
	while (my $line = <F>){
	    chomp $line;
	    my ($intwt, $intloess) = split (/\t/, $line);
	    push @{$datafin->{$k}->{$n}}, $intwt, $intloess;
	}
    }
    else {
	print STDERR "Unknown file type\n";
	die;
    }
    close (F);
}


# print
foreach my $k (sort keys %$datafin){
    
    if ($k eq "auc"){
	open (A, ">$outdir/summaries/$index/$k.coords");
	open (B, ">$outdir/summaries/$index/$k.coords_loess");
	foreach my $samp (sort {$a <=> $b} keys %{$datafin->{$k}}){
	    print A "$samp\t";
	    print B "$samp\t";
	    
	    my @aucA;
	    foreach my $nodeA (sort {$a <=> $b} keys %{$datafin->{$k}->{$samp}}){
		push (@aucA, @{$datafin->{$k}->{$samp}->{$nodeA}}[0]);
	    }
	    my $aucprintA = join "\t", @aucA;
	    print A "$aucprintA\n";
	    
	    my @aucB;
	    foreach my $nodeB (sort {$a <=> $b} keys %{$datafin->{$k}->{$samp}}){
		push (@aucB, @{$datafin->{$k}->{$samp}->{$nodeB}}[1]);
	    }
	    my $aucprintB = join "\t", @aucB;
	    print B "$aucprintB\n";
	}
	close (A);
	close (B);
    }
    elsif (($k eq "coords") or ($k eq "coords_loess")){
	open (C, ">$outdir/summaries/$index/$k.out");
	foreach my $part (sort {$a <=> $b} keys %{$datafin->{$k}}){
	    print C "$part\t";
	    
	    my @coords;
	    foreach my $nodeC (sort {$a <=> $b} keys %{$datafin->{$k}->{$part}}){
		push (@coords, $datafin->{$k}->{$part}->{$nodeC});
	    }
	    my $coordsprint = join "\t", @coords;
	    print C "$coordsprint\n";
	}
	close (C);
    }
    elsif ($k eq "integration"){
	open (I, ">$outdir/summaries/$index/$k.out");
	foreach my $nodeI (sort {$a <=> $b} keys %{$datafin->{$k}}){
	    my $intstring = join "\t", @{$datafin->{$k}->{$nodeI}};
	    print I "$nodeI\t$intstring\n";
	}
	close (I);
    }
    else {
	next;
    }
}

# print out tree annotated with fixation pts and aucs
# if there was an input tree file
if ($treefile){
    my $proj = parse (
		      '-format' => 'nexus',
		      '-file'   => $treefile,
		      '-as_project' => 1,
		      );

# NEED TO FIND RE-ROOT
#else {
#    $proj = parse (
#                   '-format' => 'newick',
#                   '-file'   => "$outdir/$index.newick.tre",
#                   '-as_project' => 1,
#                   );
#}
    my ($forest)  = @{$proj->get_forests};
    my ($tree)    = @{$forest->get_entities};
    foreach my $intnode (@{$tree->get_internals}){
	my @names;
	my @terminals = @{$intnode->get_terminals};
	for my $termnode (@terminals){
	    push (@names, $termnode->get_name);
	}
	my @srtnames = sort (@names);
	my $names = join ",", @srtnames;
	
	my $fxpt;
	if ($fixation->{$names}[1] == 1){
	    $fxpt = $fixation->{$names}[2];
	}
	else {
	    $fxpt = $fixation->{$names}[2] . "*";
	}
	
	my $auc  = $datafin->{'integration'}->{$topocode->{$names}}[0];
	
	my $annotation = $fxpt . "|" . $auc;
	$intnode->set_name($annotation);
    }
    
    open (FT, ">$outdir/summaries/$index/annotated.tre");
    my $newicktree = $tree->to_newick (nodelabels => 1);
    print FT "$newicktree\n";
    close (FT);
}

# purge if required                                                                                        
if ($purge){
    `rm -rf $outdir/cmds`;
    `rm -rf $outdir/logs`;
    `rm -rf $outdir/trees`;
    `rm -rf $outdir/output`;
}

(`rm $outdir/$index.newick.tre`) unless ($treefile); 

###SUBS###

sub round {
    my($number) = shift;
    return int($number + .5);
}

sub summarize{
    my $cfi = shift;
    my $cfitv = shift;
    my $partnum = shift;
    my $nodes = shift;
    my $type = shift;
    my $outdir = shift;
    my $topoindex = shift;
    
    # get the loess parameters
    my $span  = $cfiobj->get_span;
    my $loess = $cfiobj->get_rscript;

    # average across all samples and create data for integration
    # find starting conditions for @asym
    my @asym;
    my @x;
    my @y;
    open (AVG, ">$outdir/output/$topoindex/coords.$type");
    foreach my $include (sort {$a <=> $b} keys %$cfi){
	my $statobj = Statistics::Descriptive::Full->new();
	$statobj->add_data(@{$cfi->{$include}});
	my $mean = $statobj->mean();
	
	my $x = $include / $partnum;
	my $y = $mean / $nodes;
	push (@x, $x);
	push (@y, $y);
	
	(push @asym, $y) if ($include == 1);

	print AVG "$x\t$y\n";
    }
    close (AVG);
    
    # integrate across true data
    my $area_coords = $cfiobj->integrate_trap (\@x, \@y);
    
    # integrate the loess coordinates                                                                 
    `R --slave --args $span $outdir/output/$topoindex/coords.$type $outdir/output/$topoindex/coords_loess.$type < $loess`;
    my ($xl, $yl) = $cfiobj->parse_loess ("$outdir/output/$topoindex/coords_loess.$type");
    my $area_loess = $cfiobj->integrate_trap ($xl, $yl);

    # print the integrations
    open (INT, ">$outdir/output/$topoindex/integration.$type");
    print INT "$area_coords\t$area_loess\n";
    close (INT);
    
    # find the avg number of genes at which the asymptote is reached
    my $cacheInc;
    my $cacheX;
    my $endcondition;
    foreach my $include (sort {$b <=> $a} keys %$cfi){
	my $statobj = Statistics::Descriptive::Full->new();
        $statobj->add_data(@{$cfi->{$include}});
        my $mean = $statobj->mean();

        my $x = $include / $partnum;
        my $y = $mean / $nodes;

	if (($include - 1) == $partnum){
	    $endcondition = $y;
	    push (@asym, $endcondition);
	}
	
	if ($endcondition == 1){
	    if ($y < 1){
		push (@asym, $cacheInc, $cacheX);
		last;
	    }
	    else {
		$cacheInc = $include;
		$cacheX   = $x;
	    }
	}
	elsif ($endcondition == 0){
	    if ($y > 0){
                push (@asym, $cacheInc, $cacheX);
                last;
            }
            else {
                $cacheInc = $include;
                $cacheX   = $x;
            }
	}
	else {
	    push (@asym, "did not fix");
	    last;
	}
    }

    # if third element of @asym array is empty than that means fixation from the start:
    (push (@asym, 1, 1/$partnum)) unless ($asym[2]);
    
    # integrate each sample, true and loess
    open (AUC, ">$outdir/output/$topoindex/auc.$type");
    foreach my $sample (sort {$a <=> $b} keys %$cfitv){
	
	my @a;
	my @b;
	open (SAMP, ">$outdir/output/$topoindex/samples/$sample.$type.out");
	foreach my $include (sort {$a <=> $b} keys %{$cfitv->{$sample}}){
	    my $a = $include / $partnum;
	    my $b = $cfitv->{$sample}->{$include} / $nodes;
	    push (@a, $a);
	    push (@b, $b);
	    print SAMP "$a\t$b\n";
	}
	close (SAMP);
	
	# integrate the sample's true data
	my $area_sample = $cfiobj->integrate_trap(\@a, \@b);
	
	# integrate the sample's loess coordinates
	`R --slave --args $span $outdir/output/$topoindex/samples/$sample.$type.out $outdir/output/$topoindex/samples/$sample.$type.loess < $loess`;
	my ($xl, $yl) = $cfiobj->parse_loess ("$outdir/output/$topoindex/samples/$sample.$type.loess");
	my $area_sample_loess = $cfiobj->integrate_trap ($xl, $yl);
	
	# print the integrations                                                                             
	print AUC "$sample\t$area_sample\t$area_sample_loess\n";
    }
    close (AUC);
	  
    # print summary
    open (SUM, ">$outdir/summaries/$topoindex/summary.$type");
    foreach my $include (sort {$a <=> $b} keys %$cfi){
	my $cfistring = join "\t", @{$cfi->{$include}};
	print SUM "$include\t$cfistring\n";
    }
    close (SUM);

    return (\@asym);
}
