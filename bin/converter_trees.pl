#!/usr/bin/perl -w

# -a is the input aln file
# -f is the input format
# -b is the ouput aln file
# -o is the output format

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::TreeIO;

my ($infile, $outfile, $informat, $outformat);
GetOptions(
	   'a|infile=s'      => \$infile,
	   'b|outfile=s'     => \$outfile,
	   'f|informat=s'    => \$informat,
	   'o|outformat=s'   => \$outformat,
	   );

my $treein = new Bio::TreeIO(-file => "$infile",
			     -format => "$informat");
my $treeout = new Bio::TreeIO(-file => ">$outfile",
                                  -format => "$outformat");
while (my $tree = $treein->next_tree){
    $treeout->write_tree($tree);
}

# what follows is an intense hack necessary because the bioperl                                              
# modules wrap the nexus tree in an extra set of parens                                                      
open (TREE, "$outfile");
open (TREEO, ">$outfile.tmp");
while (my $line = <TREE>){
    chomp $line;
    if ($line =~m/tree\sBioperl\_/){
	$line =~s/\[\]//g;
	my @line = split (/\s+/, $line);
	$line[5] =~s/^.//;
	$line[5] =~s/.{2}$//;
	$line[5] .= ";";
	print TREEO "tree Bioperl_1 = [&U] $line[5]\n";
    }
    else {
	print TREEO "$line\n";
	next;
    }
}

`rm $outfile`;
`mv $outfile.tmp $outfile`;

close (TREE);
close (TREEO);
