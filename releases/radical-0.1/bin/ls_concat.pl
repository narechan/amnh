#!/usr/bin/perl -w

=head1 NAME

ls_concat.pl

=head1 SYNOPSIS

ls_concat.pl -- Calculates LS statistics.

Options:

--procs is the number of || processes
--config is the garli configuration file
    contains values for all but the following,
    which are defined as options to this program:
    datafname                                                                                               
    constraintfile                                                                                            
    ofprefix
--matrix is your data matrix
--tree is your TE treefile
--outdir is your output dir
--sets is a file with a list of charsets to analyze as a group.
    If set to 'all' you will get LS of each individual gene (not a file)

Requires the bioperl libs. 
Requires Garli to be in your path (serial version)
Requires Ls.pm

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
use Parallel::ForkManager;
use Ls;

my ($help, $procs, $configfile, $outdir, $matrix, $treefile, $sets);
GetOptions(
    'h|help'          => \$help,
    'p|procs=s'       => \$procs,
    'c|config=s'      => \$configfile,
    'o|outdir=s'      => \$outdir,
    'm|matrix=s'      => \$matrix,
    't|treefile=s'    => \$treefile,
    's|sets=s'        => \$sets,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $procs, $outdir, $matrix, $treefile, $sets){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

# make the dirs to hold the data
`mkdir -p $outdir/constraintfiles`;
`mkdir -p $outdir/data`;

#####MAIN#####

# instantiate the object and load
# default config and core matrix
my $plsobj = Ls->new;
$plsobj->load_config ($configfile);
$plsobj->load_aln    ($matrix);
$plsobj->load_tree   ($treefile);

# get the partitions, nchars, and configuration
my $charsets = $plsobj->get_charsets;
my $nchars   = $plsobj->get_nchar;
my $configuration = $plsobj->get_config;

# get the nodes and lineages
my $nodes    = $plsobj->get_leaves;
my $lineages = $plsobj->get_lineages;

# store all the taxa                                                                                      
my @alltaxa;
my $count = 0;
foreach my $node (sort keys %$nodes){
    my $newcount = @{$nodes->{$node}};

    if ($newcount > $count){
        (@alltaxa = @{$nodes->{$node}});
        $count = $newcount;
    }
}

# print node and lineage reference files
# print all possible pos and neg constraint files
my $nodeskept = {};
open (N, ">$outdir/nodes");
foreach my $node (sort keys %$nodes){

    # the constraint
    my $lineage = join (",", @{ $nodes->{$node} });
    print N "$node\t$lineage\n";

    my %node;
    foreach my $n (@{$nodes->{$node}}){
	$node{$n} = 1;
    }
    
    # the remainders
    my @remainders;
    foreach my $a (@alltaxa){
        (push (@remainders, $a)) unless (exists ($node{$a}));
    }
    my $remainders = join (",", @remainders);

    # skip the whole tree case and the root node case
    my $rcount = @remainders;
    next if ($rcount == 0);
    next if ($rcount == 1);

    $nodeskept->{$node} = 1;
    
    open (POS, ">$outdir/constraintfiles/$node.pos");
    print POS "+(($lineage),$remainders)\n";
    close (POS);
    
    open (NEG, ">$outdir/constraintfiles/$node.neg");
    print NEG "-(($lineage),$remainders)\n";
    close (NEG);
    
}
close (N);

open (L, ">$outdir/lineages");
foreach my $node (sort keys %$lineages){
    my @lins;
    foreach my $lin (sort keys %{ $lineages->{$node} }){
        push (@lins, join (",", @{ $lineages->{$node}->{$lin} }));
    }
    my $nodestream = join ("\t", @lins);
    print L "$node\t$nodestream\n";
}
close (L);

# find the group of things to analyze depending on your sets input
# if set to 'all', do the LS of individual genes
if ($sets eq "all"){

    # cycle through the partitions, and create a modified matrix,
    my $plspm = Parallel::ForkManager->new($procs);
    foreach my $charset (sort keys %$charsets){
	$plspm->start and next;
	
	# create an exclusions block file
	my ($start, $end) = split (/\-/, $charsets->{$charset});
	my $end1 = $start - 1;
	my $start2 = $end + 1;
	
	open (E, ">$outdir/data/$charset.exclusions");
	print E "Begin Assumptions;\n";
	
	if ($start == 1){
	    print E "EXSET * excluded=$start2-$nchars;\n";
	}
	elsif ($end == $nchars){
	    print E "EXSET * excluded=1-$end1;\n";
	}
	else {
	    print E "EXSET * excluded=1-$end1 $start2-$nchars;\n";
	}
	
	print E "End;\n";
	close (E);
	
	# modify the matrix with this exclusions block
	`cat $matrix $outdir/data/$charset.exclusions > $outdir/data/$charset.matrix`;
	
	# cycle through the constraints (serial here),
	opendir (DIR, "$outdir/constraintfiles");
	my @cfiles = sort (readdir (DIR));
	shift @cfiles;
	shift @cfiles;
	closedir (DIR);
	
	foreach my $cfile (@cfiles){
	    
	    # integrate matrix, constraint, and prefix changes into config
	    $configuration->{'general'}->{'datafname'} = "$outdir/data/$charset.matrix";
	    $configuration->{'general'}->{'ofprefix'} = "$outdir/data/$charset.$cfile";
	    $configuration->{'general'}->{'constraintfile'} = "$outdir/constraintfiles/$cfile";
	    
	    # output this run's specific config file
	    open (CO, ">$outdir/data/$charset.$cfile.config");
	    foreach my $header (keys %$configuration){
		print CO "[$header]\n";
		
		foreach my $key (keys %{$configuration->{$header}}){
		    print CO "$key=$configuration->{$header}->{$key}\n";
		}
	    }
	    close (CO);
	    
	    # run Garli
	    print STDERR "Running $charset $cfile\n";
	    $plsobj->run_garli("$outdir/data/$charset.$cfile.config");
	}
	
	$plspm->finish;
    }
    $plspm->wait_all_children;

    # parse the LS results                                                                                     
    open (SUM, ">$outdir/summary");
    foreach my $charset (sort keys %$charsets){
	my $coords = $charsets->{$charset};

	foreach my $node (sort {$a <=> $b} keys %$nodeskept){
	    print SUM "$charset\t$coords\t$node\t";

	    my $constlog = {};
	    foreach my $constraint ("neg", "pos"){

		open (REP, "$outdir/data/$charset.$node.$constraint.screen.log");
		while (my $line = <REP>){
		    chomp $line;

		    $line =~s/\s//g;
		    if ($line =~m/^Replicate\d\:(.*)\(best\)/){
			print SUM "$1\t";
			$constlog->{$constraint} = $1;
		    }
		    else {
			next;
		    }
		}
		close (REP);
	    }

	    my $divisor = $constlog->{"neg"};
	    ($divisor = $constlog->{"pos"}) if ($constlog->{"pos"} > $constlog->{"neg"});
	    my $ls = $constlog->{"neg"} - $constlog->{"pos"};
	    my $index = $ls / $divisor;
	    print SUM "$ls\t$index\n";
	}
    }
}

# if set to a charset list, do the LS of concat set
else {
    my $setsie = {};
    open (SETS, "$sets");
    while (my $line = <SETS>){
	chomp $line;
	$setsie->{$line} = 1;
    }
    close (SETS);

    # cycle through the charset list, and create a modified matrix,
    # excluding all those not in the file
    open (E, ">$outdir/data/exclusions");
    print E "Begin Assumptions;\n";

    my @exsets;
    foreach my $charset (sort keys %$charsets){
	if (exists ($setsie->{$charset})){ #skip the ones we want to include from ou include list!
	    next;
	}
	else {
	    push (@exsets, $charsets->{$charset});
	}
    }
    
    my $exset_string = join (" ", @exsets);
    print E "EXSET * excluded=$exset_string;\n";
    print E "End;\n";
    close (E);
	    
    # modify the matrix with this exclusions block                                                   
    `cat $matrix $outdir/data/exclusions > $outdir/data/matrix`;
	    
    # cycle through the constraints (parallel here),                                                           
    opendir (DIR, "$outdir/constraintfiles");
    my @cfiles = sort (readdir (DIR));
    shift @cfiles;
    shift @cfiles;
    closedir (DIR);

    my $plspm = Parallel::ForkManager->new($procs);
    foreach my $cfile (@cfiles){
	$plspm->start and next;

	# integrate matrix, constraint, and prefix changes into config                                       
	$configuration->{'general'}->{'datafname'} = "$outdir/data/matrix";
	$configuration->{'general'}->{'ofprefix'} = "$outdir/data/$cfile";
	$configuration->{'general'}->{'constraintfile'} = "$outdir/constraintfiles/$cfile";
	
	# output this run's specific config file                                                             
	open (CO, ">$outdir/data/$cfile.config");
	foreach my $header (keys %$configuration){
	    print CO "[$header]\n";
	    
	    foreach my $key (keys %{$configuration->{$header}}){
		print CO "$key=$configuration->{$header}->{$key}\n";
	    }
	}
	close (CO);

	# run Garli                                                                                          
	print STDERR "Running $cfile\n";
	$plsobj->run_garli("$outdir/data/$cfile.config");

	$plspm->finish;
    }
    $plspm->wait_all_children;

    # parse the LS results
    open (SUM, ">$outdir/summary");
    
    foreach my $node (sort {$a <=> $b} keys %$nodeskept){
	print SUM "$node\t";
	
	my $constlog = {};
	foreach my $constraint ("neg", "pos"){
	    
	    open (REP, "$outdir/data/$node.$constraint.screen.log");
	    while (my $line = <REP>){
		chomp $line;
		
		$line =~s/\s//g;
		if ($line =~m/^Replicate\d\:(.*)\(best\)/){
		    print SUM "$1\t";
		    $constlog->{$constraint} = $1;
		}
		else {
		    next;
		}
	    }
	    close (REP);
	}
	
	my $divisor = $constlog->{"neg"};
	($divisor = $constlog->{"pos"}) if ($constlog->{"pos"} > $constlog->{"neg"});
	my $ls = $constlog->{"neg"} - $constlog->{"pos"};
	my $index = $ls / $divisor;
	print SUM "$ls\t$index\n";
    }
}
		    
	    
	
