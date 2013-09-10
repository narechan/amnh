#!/usr/bin/perl -w

=head1 NAME

run_ild_forked.pl

=head1 SYNOPSIS

run_ild_forked.pl

Options:

--procs        is the number of procs designated
--matrix       is the data matrix
--config       is your paup config
--outdir       is your outdir for the run
--window       is your window size if doing                                                               
                LILD sildeRule (optional)                                                              
--motion       is the number of residues the                                                         
                window moves per calculation (5' -> 3')                                             
                if doing LILD slideRule (optional)   
--start        is the start position of slideRule; set to 1 by default  (optional)
--end          is the end position of slideRule; set to nchar by default    (optional)              
--file1        is your first file of partitions (optional)
--file2        is your second set of partitions (optional)
--multiple     is if your infiles have multiple charsets per ild partition (optional)

NOTE THAT THE ALL-ALL CASE AND THE SILDERULE CASE DO NOT
GENERATE THEIR OWN PAIRWISE MATRICES -- TODO!!

The matrix must be in
nexus format.

Requires the bioperl libs.
Requires Ild.pm.
Requires paup.
Requires that ild_setup.pl, ild_run.pl, and ild_parse.pl
   be in your path.

=head1 DESCRIPTION


=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
ait under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my ($help, $matrixfile, $configfile, $procs, $outdir, $window, $motion, $start, $end, $file1, $file2, $multiple);
GetOptions(
    'h|help'           => \$help,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'p|procs=s'        => \$procs,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    'a|file1=s'        => \$file1,
    'b|file2=s'        => \$file2,
    'z|multiple'     => \$multiple,
    ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $configfile, $procs, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# run the setup script
warn "Running the setup\n";
if ( ($window) and ($motion) ){
    if ( ($start) and ($end) ){
	`ild_setup.pl -m $matrixfile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end`;
    }
    else{
	`ild_setup.pl -m $matrixfile -c $configfile -o $outdir -w $window -j $motion`;
    }
}
elsif ( ($file1) and ($file2) and !$multiple){
    `ild_setup.pl -m $matrixfile -c $configfile -o $outdir -a $file1 -b $file2`;
}
elsif ( ($file1) and ($file2) and $multiple){
    `ild_setup.pl -m $matrixfile -c $configfile -o $outdir -a $file1 -b $file2 -z`;
}
elsif ( ($file1) and !($file2) ){
    `ild_setup.pl -m $matrixfile -c $configfile -o $outdir -a $file1`;
}
else {
    `ild_setup.pl -m $matrixfile -c $configfile -o $outdir`;
}

# run the run script for each experiment
opendir (E, "$outdir/cmds");
my @expts = grep (/^.+\..+$/, readdir(E));
closedir (E);

my $pm = Parallel::ForkManager->new($procs);
foreach my $expt (@expts){
    $pm->start and next;
    warn "ILD for $expt\n";
    `ild_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$expt`;
    $pm->finish;
}
$pm->wait_all_children;

# run the parse
if ( ($window) and ($motion) and ($start) and ($end) ){
    `ild_parse.pl -m $matrixfile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end`;
}
else {
    `ild_parse.pl -m $matrixfile -c $configfile -o $outdir`;
}

