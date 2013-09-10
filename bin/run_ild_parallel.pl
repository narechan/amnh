#!/usr/bin/perl -w

=head1 NAME

run_ild_parallel.pl

=head1 SYNOPSIS

run_ild_parallel.pl

Options:

--matrix       is the data matrix
--config       is your paup config
--outdir       is your outdir
--window       is your window size if doing                                                           
                LILD sildeRule (optional)                                                               
--motion       is the number of residues the                                                           
                window moves per calculation (5' -> 3')                                              
                if doing LILD slideRule (optional)                                                     
--start        is the start position of slideRule (optional)                                            
--end          is the end position of slideRule    (optional)  
--file1        is your first file of partitions (optional)                        
--file2        is your second set of partitions (optional)                        
--multiple     is if your infiles have multiple charsets per ild partition (optional)
--prefix       is your sge job name prefix

Both the matrix and the tree must be in
nexus format.

Requires the bioperl libs.
Requires Ild.pm
Requires that ild_setup.pl, ild_run.pl, and ild_parse.pl
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

my ($help, $matrixfile, $configfile, $outdir, $window, $motion, $start, $end, $prefix, $file1, $file2, $multiple);
GetOptions(
    'h|help'           => \$help,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    'w|window=s'       => \$window,
    'j|motion=s'       => \$motion,
    's|start=s'        => \$start,
    'e|end=s'          => \$end,
    'a|file1=s'        => \$file1,
    'b|file2=s'        => \$file2,
    'z|multiple'     => \$multiple,
    'p|prefix=s'       => \$prefix,	   
	   ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $configfile, $outdir, $prefix){
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
if ( ($window) and ($motion) and ($start) and ($end) ){
    `ild_setup.pl -m $matrixfile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end`;
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

# store the expts
opendir (E, "$outdir/cmds");
my @expts = grep (/^.+\..+$/, readdir(E));
closedir (E);

my $expts = @expts;

# create a shell script and submit as an SGE array job
warn "Submitting ILD array\n";
open (SHELL, ">$outdir/sge/scripts/ild.sh");
print SHELL "#!/bin/bash\n";
print SHELL "date\n";
print SHELL "source ~/.bashrc\n";
print SHELL "ild_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/\$SGE_TASK_ID.nex\n";
print SHELL "hostname\n";
print SHELL "date\n";
close (SHELL);

`qsub -r yes -N $prefix.ild -t 1-$expts:1 -S /bin/bash -cwd -o $outdir/sge/stdout/ild.out -e $outdir/sge/stderr/ild.err $outdir/sge/scripts/ild.sh`;

# run the parse when all runs have completed
open (PAR, ">$outdir/sge/scripts/parse.sh");
print PAR "#!/bin/bash\n";
print PAR "date\n";
print PAR "source ~/.bashrc\n";

if ( ($window) and ($motion) and ($start) and ($end) ){
    print PAR "ild_parse.pl -m $matrixfile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end\n";
}
else {
    print PAR "ild_parse.pl -m $matrixfile -c $configfile -o $outdir\n";
}

print PAR "hostname\n";
print PAR "date\n";
close (PAR);

warn "Submitting parse.sh\n";
`qsub -hold_jid $prefix.ild -r yes -N parse.sh -S /bin/bash -cwd -o $outdir/sge/stdout/parse.out -e $outdir/sge/stderr/parse.err $outdir/sge/scripts/parse.sh`;
