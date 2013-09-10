#!/usr/bin/perl -w

=head1 NAME

radical_parse.pl

=head1 SYNOPSIS

  radical_parse.pl -- Parse the RADICAL data.

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
my $partnum  = @charsets - 1;

# get TE topology and the root from config
my $root    = $cfiobj->get_root;

($index = "$treefile") unless ($index);
my $topoindex = {};
if ($treefile){
    $topoindex  = $cfiobj->parse_tree ($treefile, $root);
}
elsif ($dist and $index){
    open (TD, "$dist");
    while (my $line = <TD>){
	chomp $line;
	my ($num, $topo) = split (/\t/, $line);
	
	if ($num == $index){
	    my @topo = split (/\;/, $topo);
	    foreach my $node (@topo){
		my @taxa = split (/\,/, $node);
		$topoindex->{$node} = [@taxa];
	    }
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

# genrate the RADICAL stats
`mkdir -p $outdir/summaries/$index`;
`mkdir -p $outdir/output/$index/samples`;

# get trees                                                                                
opendir (L, "$outdir/trees");
my @trees = sort (readdir(L));
shift @trees;
shift @trees;
closedir (L);

# generate cfi data, and do topo comparisons
print STDERR "$index: Doing CFI and storing data\n";
my $cfi = {};
my $topo = {};
my $cfitv = {};
my $topotv = {};
foreach my $tree (@trees){
    my ($sample, $included, $nex) = split (/\./, $tree, 3);
    my $cfitopo = $cfiobj->parse_tree ("$outdir/trees/$tree", $root);
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
open (NODES, ">$outdir/summaries/$index/nodes");
foreach my $node (sort keys %$topo){
    print STDERR "$index: Sum $node\n";
    $counter++;
    print NODES "$counter\t$node\n";
    my $asymnode = 
	summarize ($topo->{$node}, $topotv->{$node}, $partnum, 1, $counter, $outdir, $index);
    $fixation->{$node} = $asymnode;
}
close (NODES);

# print out node fixation points
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

# purge if required                                                                                        
if ($purge){
    `rm -rf $outdir/cmds`;
    `rm -rf $outdir/logs`;
    `rm -rf $outdir/trees`;
    `rm -rf $outdir/output`;
}


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
