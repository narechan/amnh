#!/usr/bin/perl -w

# -m is the matrix
# -t is the tree
# -r is the root
# -l is set to 1 if you want to look at leaves as nodes
#    in addition to internal nodes; set to 0 if only internals
# -f is the input tree format
# -o is the outdir

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::TreeIO;
use TreeSupports;

my ($matrixfile, $treefile, $root, $leaves, $format, $outdir);
GetOptions(
	   'm|matrix=s'    => \$matrixfile,
	   't|tree=s'      => \$treefile,
	   'r|root=s'      => \$root,
	   'l|leaves=s'    => \$leaves,
	   'f|format=s'    => \$format,
	   'o|outdir=s'    => \$outdir,
	   );

`mkdir -p $outdir`;

#####MAIN######

# load the nodes into mem                                                                         
my $nodes = {};
my $treein = new Bio::TreeIO(-file => $treefile,
                             -format => $format);
my $treeout = new Bio::TreeIO(-file => ">$outdir/tree.index",
			      -format => $format);
my $tree = $treein->next_tree();

# re-root                                                                                             
if ($root){
    my @ns = $tree->find_node(-id => $root);
    my $n = $ns[0];
    $tree->reroot($n);
}

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
	if ($leaves == 1){
	    $nodes->{$nodeindex} = $names;
	    print N "$nodeindex\t$nodeid\n";
	}
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

# store the matrix data
my $supportobj = TreeSupports->new;
$supportobj->load_aln ($matrixfile);
my $charsets = $supportobj->get_charsets;
my $chars    = $supportobj->get_nchar;

my $alnin = Bio::AlignIO->new(-file   => $matrixfile,
                              -format => 'nexus');
my $aln = $alnin->next_aln();

my $alndata = {};
my $taxa = {};
foreach my $seq ($aln->each_seq){
    my $id        = $seq->display_id;
    $taxa->{$id} = 1;

    foreach my $charset (sort keys %$charsets){

        my $coords = $charsets->{$charset};
        my ($start, $end) = split (/\-/, $coords);

        my $partition = $seq->subseq($start, $end);
	$alndata->{$charset}->{$id} = $partition;

    }
}

# cycle through the snps and determine which are 
# diagnostic for which nodes.
# also create the all.snp report
my $genesnpcounts = {};
open (R, ">$outdir/diag.snps");
open (A, ">$outdir/all.snps");
my @taxaheader = sort keys %$taxa;
my $taxaheader = join "\t", @taxaheader;
print A "SNP\t$taxaheader\n";

foreach my $snp (keys %$alndata){
    print STDERR "$snp\n";
    $snp =~m/(\d+)/;
    my $refpos = $1;

    # print the all.snps info
    my @bases;
    foreach my $taxon (sort keys %{$alndata->{$snp}}){
	push (@bases, $alndata->{$snp}->{$taxon});
    }
    my $bases = join "\t", @bases;
    print A "$snp\t$bases\n";
    
    # cycle through nodes
    foreach my $nodeindex (sort {$a <=> $b} keys %$nodes){
	my @nodein;
	my @nodeout;
	
	# populate the outgroup
	my @nodeouttaxa;
	foreach my $taxon (keys %$taxa){
	    if (exists ($nodes->{$nodeindex}->{$taxon})){
		next;
	    }
	    else {
		push (@nodeout, $alndata->{$snp}->{$taxon});
		push (@nodeouttaxa, $taxon);
	    }
	}
	my $nodeoutstring = join ",", sort (@nodeouttaxa);
       
	
	# populate the ingroup
	foreach my $leaf (keys %{$nodes->{$nodeindex}}){
	    push (@nodein, $alndata->{$snp}->{$leaf});
	}
	my $nodeinstring = join ",", sort (keys %{$nodes->{$nodeindex}});

	# check to see if the ingroup is internally consistent
	if (keys %{{map {$_, 1} @nodein}} == 1){
	    my $basein = $nodein[0];
	    
	    # check to see if the outgroup bases are all diff
	    # from the ingroup -- don't have to be internally consistent
	    my $compres = 0;
	    my $baseouts = {};
	    foreach my $baseout (@nodeout){
		if ($baseout eq $basein){
		    $compres++;
		}
		else{
		    $baseouts->{$baseout}++;
		    next;
		}
	    }
	    
	    # print if the conditions are satisfied
	    if ($compres > 0){
		next;
	    }
	    else {
		my $baseoutstr = join (",", sort (keys %$baseouts));
		print R "$snp\t$nodeindex\t$basein\t$baseoutstr\n";
	    }
	}
	else {
	    next;
	}
    }
}
close (R);
close (A);
