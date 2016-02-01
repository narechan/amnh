#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Radical;

#my ($help, $treefile, $listfile, $root);
my ($help, $treefile, $listfile, $rootfile);
GetOptions(
    'h|help'          => \$help,
    't|treefile=s'   => \$treefile,
    'l|listfile=s'       => \$listfile,
#    'r|root=s' => \$root,
    'r|rootfile=s'   => \$rootfile,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $listfile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

print STDERR "Working on $treefile\n";

# store the list
my $list = {};
open (L, "$listfile");
while (my $line = <L>){
    chomp $line;
    $list->{$line} = 1;
}
close (L);

# store the roots from the file
my $roots = {};
open (R, "$rootfile");
while (my $line = <R>){
    chomp $line;
    my ($fam, $sp, $seq) = split (/\t/, $line);
    my $root = $sp . "#" . $seq;
    $roots->{$fam} = $root;
}
close (R);

# get the treefile's OID group
my @treefile = split (/\./, $treefile);

# instantiate the object and load stuff we need
my $cfiobj = Radical->new;

# get the tree topology
#my $topo = $cfiobj->parse_tree ($treefile, $root, 'newick');
#my $topo = $cfiobj->parse_tree ($treefile, $roots->{$treefile[3]}, 'newick');
my $topo = $cfiobj->parse_tree ($treefile, $roots->{$treefile[1]}, 'newick');

my $seqskept  = {};
my $nodeskept = {};
foreach my $node (keys %$topo){
    my @seqs = split (/\,/, $node);
    my $taxa = {};
    my $counter = 0;
    
    foreach my $seq (@seqs){
	my ($taxon, $acc) = split (/\#/, $seq);
	$taxa->{$taxon}++;
	$counter++;
    }
    
    my $signal = 0;
    foreach my $t (keys %$taxa){
	if (exists ($list->{$t})){
	    next;
	}
	else {
	    $signal = 1;
	    last;
	}
    }
    
    if ($signal == 0){
	$nodeskept->{$node} = $counter;
	
	my @seqskept = split (/\,/, $node);
	foreach my $seqk (@seqskept){
	    $seqskept->{$seqk} = 1;
	}
    }

    else {
	next;
    }
}

if (%$nodeskept){
    my $cladeskept = {};

    foreach my $seqkept (keys %$seqskept){
	my $topnode;
	my $topnodecount = 0;

	foreach my $node (keys %$nodeskept){
	    my $signal = 0;

	    my @s = split (/\,/, $node);
	    foreach my $s (@s){
		if ($s eq $seqkept){
		    $signal = 1;
		}
		else {
		    next;
		}
	    }
	    
	    if ($signal == 1){
		if ($nodeskept->{$node} > $topnodecount){
		    $topnodecount = $nodeskept->{$node};
		    $topnode      = $node;
		}
		else {
		    next;
		}
	    }
	    else {
		next;
	    }
	}
	$cladeskept->{$topnode} = 1;
	
    }
    my $cladecounter = 0;
    foreach my $clade (keys %$cladeskept){
	$cladecounter++;
	print "$treefile\tclade$cladecounter\t$clade\n";
    }
}
else {
    print "$treefile\tNONE\n";
}
