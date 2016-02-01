#!/usr/bin/perl -w

=head1 NAME

run_lild_parallel.pl

=head1 SYNOPSIS

run_lild_parallel.pl

Options:

--matrix       is the data matrix
--tree         is the total evidence tree    
--config       is your paup config
--outdir       is your outdir
--window       is your window size if doing                                                            
                LILD sildeRule (optional)                                                            
--motion       is the number of residues the                                                         
                window moves per calculation (5' -> 3')                                              
                if doing LILD slideRule (optional)
--start        is the start position of slideRule (optional)                      
--end          is the end position of slideRule    (optional)

Both the matrix and the tree must be in
nexus format.

Requires the bioperl libs.
Requires Lild.pm
Requires that lild_setup.pl, lild_run.pl, and lild_parse.pl
   be in your path.

=head1 DESCRIPTION

This script runs in batch mode without the use of SGE
job arrays.

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $treefile, $matrixfile, $configfile, $outdir, $window, $motion, $start, $end);
GetOptions(
    'h|help'           => \$help,
    't|tree=s'         => \$treefile,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $matrixfile, $configfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

# make sge specific dirs
`mkdir -p $outdir/sge/scripts`;
`mkdir -p $outdir/sge/stdout`;
`mkdir -p $outdir/sge/stderr`;

#####MAIN#####

# run the setup script
warn "Running the setup\n";
if ( ($window) and ($motion) ){
    `lild_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end`;
}
else {
    `lild_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir`;
}

# store the expts
opendir (E, "$outdir/expts");
my @expts = grep (/^.+\..+$/, readdir(E));
closedir (E);

# create a shell script for every expt and submit
foreach my $expt (@expts){
    warn "Submitting LILD for lild.$expt\n";
    
    open (SHELL, ">$outdir/sge/scripts/$expt.sh");
    print SHELL "#!/bin/bash\n";
    print SHELL "date\n";
    print SHELL "source ~/.bashrc\n";
    print SHELL "lild_run.pl -m $matrixfile -c $configfile -o $outdir -t $treefile -e $outdir/expts/$expt\n";
    print SHELL "hostname\n";
    print SHELL "date\n";
    close (SHELL);

    `qsub -r yes -N lild.$expt -S /bin/bash -cwd -V -o $outdir/sge/stdout/$expt.out -e $outdir/sge/stderr/$expt.err $outdir/sge/scripts/$expt.sh`;
} 

# run the parse when all runs have completed
open (PAR, ">$outdir/sge/scripts/parse.sh");
print PAR "#!/bin/bash\n";
print PAR "date\n";
print PAR "source ~/.bashrc\n";
print PAR "lild_parse.pl -m $matrixfile -c $configfile -o $outdir -t $treefile\n";
print PAR "hostname\n";
print PAR "date\n";
close (PAR);

warn "Submitting parse.sh\n";
`qsub -hold_jid lild.* -r yes -N parse.sh -S /bin/bash -cwd -o $outdir/sge/stdout/parse.out -e $outdir/sge/stderr/parse.err $outdir/sge/scripts/parse.sh`;
