#!/usr/bin/perl -w

#####SETUP#####
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::TreeIO;
use Bio::Tree::RandomFactory;
use Algorithm::Numerical::Shuffle qw /shuffle/;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $taxanum, $treenum, $genenum, $outdir, $genelen, $treesamp, $genesamp, $fracmissing);
# treenum and genenum correspond to the number of trees and genes to build
# treesamp and genesamp correspond to the number of trees and genes to sample
# fracmissing should translate into whole numbers of missing taxa, no fractions
GetOptions(
    'h|help'          => \$help,
    't|taxa=s'        => \$taxanum,
    'x|trees=s'       => \$treenum,
    'g|genes=s'       => \$genenum,
    'o|outdir=s'      => \$outdir,
    'l|genelen=s'      => \$genelen,
    'y|genesamp=s'    => \$genesamp,
    'z|treesamp=s'    => \$treesamp,
    'm|missing=s'     => \$fracmissing,
	   ) or pod2usage;
pod2usage if $help;

`mkdir -p $outdir/fastas`;

#####MAIN#####

# generate the taxon names
my @taxa;
for (my $j = 0; $j < $taxanum; $j++){
    my $taxon = "sp" . $j;
    push (@taxa, $taxon);
}

# initialize a TreeIO writer
my $out = Bio::TreeIO->new(-format => "newick",
                           -file   => ">$outdir/randomtrees.tre");
my $factory = new Bio::Tree::RandomFactory(-taxa => \@taxa);
 
# generate the random trees
for (my $i = 0; $i < $treenum; $i++) {
    $out->write_tree($factory->next_tree);
}

# parse through the trees and store
my $treefile = "$outdir/randomtrees.tre";
open my $th, '<', $treefile;
$/ = undef;
my $treestring = <$th>;
close $th;
$treestring =~s/Node\d+//g;
my @trees = split (/\;/, $treestring);
`rm $outdir/randomtrees.tre`;

# build the number of genes requested for each tree
my $data = {};
my $tcounter = -1;
foreach my $tree (@trees){
    $tcounter++;
    my $prand = int(rand(1000000000000000000));
    
    # print the tree
    open (T, ">$outdir/$tcounter.tre");
    print T "$tree\n";
    close (T);

    # run seq-gen for the tree to generate the 
    # max number of genes that can be sampled from this tree
    `seq-gen -z$prand -mJTT -n$genenum -l$genelen < $outdir/$tcounter.tre > $outdir/$tcounter.data`;
   
    # parse the alignments
    open my $ah, '<', "$outdir/$tcounter.data";
    my $alnstring = <$ah>;
    close $ah;
    
    # print the alignments
    my $acounter = -1;
    my @alns = split (/\n\s/, $alnstring);
    foreach my $aln (@alns){
	$acounter++;
	open (A, ">$outdir/$acounter.$tcounter.phy");
	(chop $aln) if ($aln =~m/\n$/);
	print A "$aln\n";
	close (A);
    
	# transform the alignments
	open (B, ">$outdir/$acounter.$tcounter.fa");
	open my $xh, '<', "$outdir/$acounter.$tcounter.phy";
	my $phystring = <$xh>;
	close $xh;
	my @phystring = split (/\n/, $phystring);
	shift @phystring;
	foreach my $line (@phystring){
	    my ($sp, $seq) = split (/\s+/, $line);
	    print B ">$sp\n$seq\n";
	}
	close (B);
	
	# suck the alignment into a datastruc
	open my $fh, '<', "$outdir/$acounter.$tcounter.fa";
        my $fastring = <$fh>;
        close $fh;
	$data->{$tcounter}->{$acounter} = $fastring;
    }
}

# purge          
`rm $outdir/*.tre`;
`rm $outdir/*.data`;
`rm $outdir/*.phy`;
`rm $outdir/*.fa`;

# ramdomly pick the trees requested
my @treekeys = keys %$data;
my $treecache = {};
my @treeskept;
for (my $t = 1; $t <= $treesamp; $t++){
    my $treeindex = int (rand @treekeys);
    if (exists ($treecache->{$treeindex})){
	($treeindex = int (rand @treekeys)) until (! exists($treecache->{$treeindex}));
    }
    $treecache->{$treeindex} = 1;
    push (@treeskept, $treeindex);
}

# randomly pick the genes among the trees selected
my $finalpicks = {};
my $genecache = {};
for (my $g = 1; $g <= $genesamp; $g++){
    my $randtree = $treeskept[int (rand @treeskept)];
    my @genekeys = keys %{$data->{$randtree}};
    my $geneindex = int (rand @genekeys);
    if (exists($genecache->{$randtree}->{$geneindex})){
	($geneindex = int (rand @genekeys)) until (! exists($genecache->{$randtree}->{$geneindex}));
    }
    $genecache->{$randtree}->{$geneindex} = 1;
    $finalpicks->{$randtree}->{$geneindex} = $data->{$randtree}->{$geneindex};
}

# print the picks
my $missingTaxa = $taxanum * $fracmissing; # want this to be a whole number otherwise this hack with break
foreach my $a (keys %$finalpicks){
    foreach my $b (keys %{$finalpicks->{$a}}){
	my $filename = $a . "_" . $b;
	open (WH, ">$outdir/fastas/$filename");
	print WH "$finalpicks->{$a}->{$b}";
	close (WH);
	
	my @shuffledtaxa = shuffle @taxa;
	my @missings = splice @shuffledtaxa, 0, $missingTaxa;
	my $missings = {};
	foreach my $missing (@missings){
	    $missings->{$missing} = 1;
	}
	
	open (M, ">$outdir/fastas/$filename.missings");
	my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$outdir/fastas/$filename");
	while (my $sequence_obj = $seqin->next_seq()){
	    my $id       = $sequence_obj->display_id();
	    my $seq      = $sequence_obj->seq();
	    my $len      = length ($seq);
	    
	    if (exists ($missings->{$id})){
		print M ">$id\n";
		print M '?' x $len;
		print M "\n";
	    }
	    else{
		print M ">$id\n$seq\n";
	    }
	}
	close (M);
	`rm $outdir/fastas/$filename`;
	`mv $outdir/fastas/$filename.missings $outdir/fastas/$filename`;
    }
}


    
=head
    # parse the resulting genes into separate fasta alignments
    my $alnin = Bio::AlignIO->new(-file   => "$outdir/$tcounter.data",
                                  -format => "phylip");
    my $acounter = 0;
    while (my $alnobj = $alnin->next_aln){ 
	$acounter++;
	open (A, ">$outdir/$acounter.$tcounter.aln");
	foreach my $seq ($alnobj->each_seq){
	    my $id        = $seq->display_id;
	    my $sequence = $seq->seq;
	    print A ">$id\n$sequence\n";
	}
    }
=cut


    
    

=head
# build the partitions for seq-gen
my $titerator = $genenum / $treenum;
my $tindex = 0;
my $counter = -1;
my $gcounter = -1;
until ($gcounter == ($genenum - 1)){
    $gcounter++;
    $counter++;
    if ($counter == $titerator){
	$tindex++;
	$counter = 0;
	print STDERR "Loop1\t$gcounter\t$counter\t$tindex\n";
	print "[$genelen]$trees[$tindex]\n";
    }
    else {
	print STDERR "Loop2\t$gcounter\t$counter\t$tindex\n";
	print "[$genelen]$trees[$tindex]\n";
    }
}
=cut
