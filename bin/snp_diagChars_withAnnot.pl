#!/usr/bin/perl -w

# -m is the matrix
# -t is the tree
# -r is the root
# -l is set to 1 if you want to look at leaves as nodes
#    in addition to internal nodes; set to 0 if only internals
# -f is the input tree format
# -o is the outdir
# -g is the gff3 annotation file
# -a is the tuberculist annotation file
# -x is the reference sequence id

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::TreeIO;
use TreeSupports;

my ($matrixfile, $treefile, $root, $leaves, $format, $outdir, $gff3file, $annotfile, $refid);
GetOptions(
	   'm|matrix=s'    => \$matrixfile,
	   't|tree=s'      => \$treefile,
	   'r|root=s'      => \$root,
	   'l|leaves=s'    => \$leaves,
	   'f|format=s'    => \$format,
	   'o|outdir=s'    => \$outdir,
	   'g|gff3=s'      => \$gff3file,
	   'a|annotfile=s' => \$annotfile,
	   'x|refid=s'     => \$refid,
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

# store tuberculist
my $listdata = {};
open (F, "$annotfile");
while (my $line = <F>){
    chomp $line;
    my @data = split (/\t/, $line);
    $listdata->{$data[1]} = $data[0];
}
close (F);

# store gff3
my $gff3str = {};
my $gff3att = {};
my $gff3typ = {};

open (G, "$gff3file");
while (my $line = <G>){
    next if ($line =~m/^\#/);
    chomp $line;
    
    my ($chrom,
	$source,
	$type,
	$start,
	$end,
	$score,
	$strand,
	$phase,
	$attr) = split (/\t/, $line);
    
    my @attrs = split (/\;/, $attr);
    
    # store the information for tRNA, rRNA, CDS (mRNA),                       
    # and transcript only                                                     
    if ($type =~m/CDS|tRNA|rRNA|transcript/){
	my $range = $start . "-" . $end;
	$gff3str->{$chrom}->{$range} = $strand;
	$gff3att->{$chrom}->{$range} = [@attrs];
	$gff3typ->{$chrom}->{$range} = $type;
    }
    else {
	next;
    }
}
close (G);

# cycle through the snps and determine which are 
# diagnostic for which nodes.
# also create the all.snp report
my $genesnpcounts = {};
open (R, ">$outdir/diag.snps");
open (A, ">$outdir/all.snps");
my @taxaheader = sort keys %$taxa;
my $taxaheader = join "\t", @taxaheader;
print A "SNP\tGENE\tSNPCOUNT\t$taxaheader\n";

foreach my $snp (keys %$alndata){
    $snp =~m/(\d+)/;
    my $refpos = $1;

    # determine which gene this snp sits in
    my $genename;
    my $genedesc;
    foreach my $range (sort keys %{$gff3str->{$refid}}){
        my ($start, $end) = split (/\-/, $range);

	# Get the gene annotations -- selected vars are empty if              
	# key/val pairs don't exist in the attr                               
	my @attrs = @{$gff3att->{$refid}->{$range}};
	my $name;
	my $desc;
	foreach my $attr (@attrs){
	    my ($key, $value) = split (/\=/, $attr);
	    if ($key eq "Name"){
		$name = $value;
	    }
	    if ($key eq "description"){
		$desc = $value;
	    }
	}

	
        # if it's in a gene, find out which one, pull the sequence,               
        # and annotate as syn or nonsyn                                           
        if (($start <= $refpos) and ($refpos <= $end)){
	    $genesnpcounts->{$name}++;
	    $genename = $name;
	    $genedesc = $desc;
	    last;
	}
    }
    
    # print the all.snps info
    my @bases;
    my $snpcount = 0;
    foreach my $taxon (sort keys %{$alndata->{$snp}}){
	($snpcount++) if ($alndata->{$snp}->{$refid} ne $alndata->{$snp}->{$taxon});
	push (@bases, $alndata->{$snp}->{$taxon});
    }
    my $bases = join "\t", @bases;
    if ($genename){
	print A "$snp\t$genename\t$snpcount\t$bases\n";
    }
    else {
	print A "$snp\tINTERGENIC\t$snpcount\t$bases\n";
    }
    
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
		
		if ($genename){
		    print R "$snp\t$nodeindex\t$basein\t$baseoutstr\t$genename\t";
		    
		    if ($genedesc){
			print R "$genedesc\t";
		    }
		    else {
			print R "NO DESC\t";
		    }
		    
		    if (exists ($listdata->{$genename})){
			print R "$listdata->{$genename}\n";
		    }
		    else {
			print R "NO FUNC CAT\n";
		    }
		}
		else {
		    print R "$snp\t$nodeindex\t$basein\t$baseoutstr\tINTERGENIC\n";
		}
	    }
	}
	else {
	    next;
	}
    }
}
close (R);
close (A);

# print out snps per gene
open (G, ">$outdir/gene.snps");
foreach my $range (sort keys %{$gff3str->{$refid}}){
    my ($start, $end) = split (/\-/, $range);
    
    my @attrs = @{$gff3att->{$refid}->{$range}};
    my $name;
    my $desc;
    foreach my $attr (@attrs){
	my ($key, $value) = split (/\=/, $attr);
	if ($key eq "Name"){
	    $name = $value;
	}
	if ($key eq "description"){
	    $desc = $value;
	}
    }

    print G "$name\t$start\t$end\t";
    if ($desc){
	print G "$desc\t";
    }
    else {
	print G "NO DESC\t";
    }
    
    if (exists ($listdata->{$name})){
	print G "$listdata->{$name}\t";
    }
    else {
	print G "NO FUNC CAT\t";
    }
    
    if (exists ($genesnpcounts->{$name})){
	print G "$genesnpcounts->{$name}\n";
    }
    else {
	print G "NO SNPS\n";
    }
}
close (G);

# print out snp annotations
