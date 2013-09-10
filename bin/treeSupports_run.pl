#!/usr/bin/perl -w

=head1 NAME

treeSupports_run.pl

=head1 SYNOPSIS

treeSupports_run.pl

Options:

--experiment is your charset/node pair file
--support is your experiment type (bs or pbs)    
--config is your run configuration
--matrix is your alignfile
--tree is your treefile
--kind is the kind of bs node test you want to run
    1 is one informative taxa per lineage
    0 is two informative taxa per node
--outdir is your outdir for the run
--window       is your window size if doing                                                                
    LILD sildeRule (optional)                                                                         
--motion       is the number of residues the                                                          
    window moves per calculation (5' -> 3')                                                     
    if doing LILD slideRule (optional) 
--start        is the start position of slideRule                                                        
    (optional)                                                                                         
--end          is the end position of slideRule                                                        
    (optional)                             

Both the matrix and tree must be in 
nexus format.

Dependencies:

PAUP must be installed and in $PATH
Requires the bioperl libs. 
Requires TreeSupports.pm

=head1 DESCRIPTION


=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2012 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use lib '/Users/apurvawork/amnh/bin/treeSupports';

use strict;
use Getopt::Long;
use Pod::Usage;
use TreeSupports;

my ($help, $exptfile, $support, $matrixfile, $treefile, $configfile, $outdir, $kind, $window, $motion, $start, $end);
GetOptions(
    'h|help'          => \$help,
    'f|experiment=s'  => \$exptfile,
    'x|support=s'     => \$support,
    'm|matrix=s'      => \$matrixfile,
    't|treefile=s'    => \$treefile,
    'c|config=s'      => \$configfile,
    'o|outdir=s'      => \$outdir,
    'k|kind=s'        => \$kind,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    ) or pod2usage;
pod2usage if $help;

for my $option ($exptfile, $support, $matrixfile, $treefile, $configfile, $outdir, $kind){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $supportobj = TreeSupports->new;
$supportobj->load_config ($configfile);
my $root = $supportobj->get_root;
my $config = $supportobj->get_config;

$supportobj->load_aln    ($matrixfile);
$supportobj->load_tree   ($treefile, $root, 'nexus');

# override partitions in the nxs file if slideRule chosen
if ( ($window) and ($motion) and ($start) and ($end) ){
    $supportobj->generate_partitions ($window, $motion, $start, $end);
}

# create and execute the paup command file for bs
if ($support eq "ndi"){
    my $node;
    my $constraint;
    open (EXPT, "$exptfile");
    while (my $line = <EXPT>){
        chomp $line;
        ($node, $constraint) = split (/\t/, $line);
    }
    close (EXPT);

    my @array = ("without", "with");
    foreach my $enforce (@array){
	my $noderun =
	    $supportobj->build_ndi (
				    $matrixfile,
				    $node,
				    "$outdir/ndi/cmds",
				    "$outdir/ndi/logs",
				    "$outdir/ndi/trees",
				    $enforce,
				    );
	
	$supportobj->run_paup (
			       "$outdir/ndi/cmds",
                               "ndi",
			       $node,
			       $enforce,
                               );
    }
}
    
elsif ($support eq "bs"){
    my $charset;
    my $coords;
    open (EXPT, "$exptfile");
    while (my $line = <EXPT>){
	chomp $line;
	($charset, $coords) = split (/\t/, $line); 
    }
    close (EXPT);

    my $full_taxa = {};
    my $all_taxa  = {};
    my $missing_taxa = {};
    
    # a way to short circuit the bs missing taxa check if taxa are defined in the config
    if ((exists ($config->{'ALLTAXA'}))){
	my @fulltaxa = split (/\s/, $config->{'FULLTAXA'});
	my @alltaxa  = split (/\s/, $config->{'ALLTAXA'});
	my @missingtaxa = split (/\s/, $config->{'MISSINGTAXA'});
	
	foreach my $ftaxa (@fulltaxa){
	    $full_taxa->{$ftaxa} = 1;
	}
	foreach my $ataxa (@alltaxa){
	    $all_taxa->{$ataxa} = 1;
	}
	foreach my $mtaxa (@missingtaxa){
	    $missing_taxa->{$mtaxa} = 1;
	}
    }
    
    else {
	$supportobj->transform_aln (
				    $matrixfile, 
				    "nexus", 
				    "fasta", 
				    "$outdir/bs/cmds", 
				    $charset, 
				    $coords
				    );
	
	($full_taxa, $all_taxa, $missing_taxa) = 
	    $supportobj->find_inf_taxa (
					$matrixfile, 
					"fasta", 
					"$outdir/bs/cmds", 
					$charset
					);
    }
    
    # check to see if the topology of the 
    # PBS without tree is the same as the topology of 
    # the SA tree with respect to informative taxa.
    # Store the information.
    
    my $partition = 
	$supportobj->check_topo_build_bs (
	    $matrixfile,
	    $treefile,
	    $full_taxa, 
	    $all_taxa, 
	    $missing_taxa,			  
	    "$outdir/bs/cmds", 
	    "$outdir/bs/logs",
	    "$outdir/bs/prunes",
	    "$outdir/pbs/trees",			  
	    $charset,
            $coords, 
	    $kind					  
	);
    
    $supportobj->run_paup (
	"$outdir/bs/cmds", 
	"bs", 
	$partition
	); 
}

elsif ($support eq "pbs"){
    my $node;
    my $constraint;
    open (EXPT, "$exptfile");
    while (my $line = <EXPT>){
	chomp $line;
	($node, $constraint) = split (/\t/, $line);
    }
    close (EXPT);

    my @array;
    if ($node eq "unconstrained"){
	@array = ("none");
    }
    else {
	@array = ("without", "with");
    }

    foreach my $enforce (@array){
	my $noderun = 
	    $supportobj->build_pbs (
				    $matrixfile,
				    $node,
				    "$outdir/pbs/cmds",
				    "$outdir/pbs/logs",
				    "$outdir/pbs/trees",
				    $enforce
				    );
	
	$supportobj->run_paup (
			       "$outdir/pbs/cmds",
			       "pbs",
			       $node,
			       $enforce
			       );
    }
}

else{
    print STDERR "Expt type not recognized\n";
    die;
}
