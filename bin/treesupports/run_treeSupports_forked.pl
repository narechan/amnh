#!/usr/bin/perl -w

=head1 NAME

run_treeSupports_forked.pl

=head1 SYNOPSIS

run_treeSupports_forked.pl

Options:

--procs        is the max number of procs to fork
--matrix       is the data matrix
--tree         is the total evidence tree    
--config       is your paup config, and your experiments config
--outdir       is your outdir
--kind         is the kind of bs node test you want to run
    1 is one informative taxa per lineage                                                                
    2 is two informative taxa per node 
--window       is your window size if doing                                                            
    sildeRule (optional)                                                                      
--motion       is the number of residues the                                                           
    window moves per calculation (5' -> 3')                                                      
    if doing slideRule (optional)
--start        is the start position of slideRule                                                          
    (optional)                                                                                           
--end          is the end position of slideRule                                                   
    (optional)                             

Both the matrix and the tree must be in
nexus format.

Dependencies:

PAUP must be installed and in $PATH
Requires the bioperl libs.
Requires TreeSupports.pm

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

use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my ($help, $procs, $treefile, $matrixfile, $configfile, $outdir, $kind, $window, $motion, $start, $end);
GetOptions(
    'h|help'           => \$help,
    'p|procs=s'        => \$procs,   
    't|tree=s'         => \$treefile,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    'k|kind=s'         => \$kind,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    ) or pod2usage;
pod2usage if $help;

for my $option ($procs, $treefile, $matrixfile, $configfile, $outdir, $kind){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# run the setup script
warn "Running the setup\n";
if ( ($window) and ($motion) and ($start) and ($end) ){
    `treeSupports_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end`;
}
else {
    `treeSupports_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir`;
}

# run the run script for each NDI experiment                                                             
opendir (NDI, "$outdir/ndi/expts");
my @ndiexpts = grep (/^.+\..+$/, readdir(NDI));
closedir (NDI);

my $pmndi = Parallel::ForkManager->new($procs);
foreach my $ndiexpt (@ndiexpts){
    $pmndi->start and next;
    warn "NDI for $ndiexpt\n";
    if ( ($window) and ($motion) and ($start) and ($end) ){
        `treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x ndi -k $kind -o $outdir -f $outdir/ndi/expts/$ndiexpt -w $window -j $motion -s $start -e $end`;
    }
    else{
        `treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x ndi -k $kind -o $outdir -f $outdir/ndi/expts/$ndiexpt`;
    }
    $pmndi->finish;
}
$pmndi->wait_all_children;


# run the run script for each PBS experiment                                                   
opendir (PBS, "$outdir/pbs/expts");
my @pbsexpts = grep (/^.+\..+$/, readdir(PBS));
closedir (PBS);

my $pmpbs = Parallel::ForkManager->new($procs);
foreach my $pbsexpt (@pbsexpts){
    $pmpbs->start and next;
    warn "PBS for $pbsexpt\n";
    if ( ($window) and ($motion) and ($start) and ($end) ){
	`treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x pbs -k $kind -o $outdir -f $outdir/pbs/expts/$pbsexpt -w $window -j $motion -s $start -e $end`;
    }
    else{
	`treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x pbs -k $kind -o $outdir -f $outdir/pbs/expts/$pbsexpt`;
    }
    $pmpbs->finish;
}
$pmpbs->wait_all_children;

# run the script foreach BS expt
opendir (BS, "$outdir/bs/expts");
my @bsexpts = grep (/^.+\..+$/, readdir(BS));
closedir (BS);

my $pmbs = Parallel::ForkManager->new($procs);
foreach my $bsexpt (@bsexpts){
    $pmbs->start and next;
    warn "BS for $bsexpt\n";
    if ( ($window) and ($motion) and ($start) and ($end) ){
	`treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x bs -k $kind -o $outdir -f $outdir/bs/expts/$bsexpt -w $window -j $motion -s $start -e $end`;
    }
    else{
	`treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x bs -k $kind -o $outdir -f $outdir/bs/expts/$bsexpt`;
    }
    $pmbs->finish;
} 
$pmbs->wait_all_children;

# run the parse
if ( ($window) and ($motion) and ($start) and ($end) ){
    `treeSupports_parse.pl -m $matrixfile -t $treefile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end`;
}
else{
    `treeSupports_parse.pl -m $matrixfile -t $treefile -c $configfile -o $outdir`;
}
