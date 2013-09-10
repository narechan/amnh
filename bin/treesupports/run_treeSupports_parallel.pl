#!/usr/bin/perl -w

=head1 NAME

run_treeSupports_parallel.pl

=head1 SYNOPSIS

run_treeSupports_parallel.pl

Options:

--matrix       is the data matrix
--tree         is the total evidence tree    
--config       is your paup config
--outdir       is your outdir
--kind         is the kind of bs node test you want to run                        
    1 is one informative taxa per lineage                                  
    2 is two informative taxa per node                                         
--window       is your window size if doing                                 
    LILD sildeRule (optional)                                               
--motion       is the number of residues the                                
    window moves per calculation (5' -> 3')                                      
    if doing LILD slideRule (optional)                                            
--start        is the start position of slideRule                                
    (optional)                                                                   
--end          is the end position of slideRule                                  
    (optional)                                       

Both the matrix and the tree must be in
nexus format.

Requires the bioperl libs.
Requires TreeSupports.pm
Requires that treeSupports_setup.pl, treeSupports_run.pl, and treeSupports_parse.pl
   be in your path.

=head1 DESCRIPTION


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

my ($help, $treefile, $matrixfile, $configfile, $outdir, $kind, $window, $motion, $start, $end);
GetOptions(
    'h|help'           => \$help,
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

for my $option ($treefile, $matrixfile, $configfile, $outdir, $kind){
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
    `treeSupports_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end`;
}
else {
    `treeSupports_setup.pl -t $treefile -m $matrixfile -c $configfile -o $outdir`;
}

# run the run script for each PBS experiment                                                           
opendir (PBS, "$outdir/pbs/expts");
my @pbsexpts = grep (/^.+\..+$/, readdir(PBS));
closedir (PBS);

foreach my $pbsexpt (@pbsexpts){
    warn "PBS for $pbsexpt\n";

    open (SHELL, ">$outdir/sge/scripts/$pbsexpt.pbs.sh");
    print SHELL "#!/bin/bash\n";
    print SHELL "date\n";
    print SHELL "source ~/.bashrc\n";
    if ( ($window) and ($motion) and ($start) and ($end) ){
        print SHELL "treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x pbs -k $kind -o $outdir -f $outdir/pbs/expts/$pbsexpt -w $window -j $motion -s $start -e $end\n";
    }
    else{
        print SHELL "treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x pbs -k $kind -o $outdir -f $outdir/pbs/expts/$pbsexpt\n";
    }
    print SHELL "hostname\n";
    print SHELL "date\n";
    close (SHELL);
    
    `qsub -r yes -N ts.pbs.$pbsexpt -S /bin/bash -cwd -o $outdir/sge/stdout/$pbsexpt.out -e $outdir/sge/stderr/$pbsexpt.err $outdir/sge/scripts/$pbsexpt.pbs.sh`;
}


# run the run for each BS expt
# BS expts rely on trees from pbs -- so all these must run only
# after PBS runs are complete
opendir (BS, "$outdir/bs/expts");
my @bsexpts = grep (/^.+\..+$/, readdir(BS));
closedir (BS);

# create shell scripts for each BS expt and submit
foreach my $bsexpt (@bsexpts){
    warn "BS for $bsexpt\n";
    
    open (SHELL, ">$outdir/sge/scripts/$bsexpt.bs.sh");
    print SHELL "#!/bin/bash\n";
    print SHELL "date\n";
    print SHELL "source ~/.bashrc\n";
    if ( ($window) and ($motion) and ($start) and ($end) ){
	print SHELL "treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x bs -k $kind -o $outdir -f $outdir/bs/expts/$bsexpt -w $window -j $motion -s $start -e $end\n";
    }
    else{
        print SHELL "treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x bs -k $kind -o $outdir -f $outdir/bs/expts/$bsexpt\n";
    }
    print SHELL "hostname\n";
    print SHELL "date\n";
    close (SHELL);

    `qsub -hold_jid ts.pbs.* -r yes -N ts.bs.$bsexpt -S /bin/bash -cwd -o $outdir/sge/stdout/$bsexpt.out -e $outdir/sge/stderr/$bsexpt.err $outdir/sge/scripts/$bsexpt.bs.sh`;
}


# run the parse
open (PAR, ">$outdir/sge/scripts/parse.sh");
print PAR "#!/bin/bash\n";
print PAR "date\n";
print PAR "source ~/.bashrc\n";
if ( ($window) and ($motion) and ($start) and ($end) ){
    print PAR "treeSupports_parse.pl -m $matrixfile -t $treefile -c $configfile -o $outdir -w $window -j $motion -s $start -e $end\n";
}
else{
    print PAR "treeSupports_parse.pl -m $matrixfile -t $treefile -c $configfile -o $outdir\n";
}
print PAR "hostname\n";
print PAR "date\n";
close (PAR);

`qsub -hold_jid ts.* -r yes -N parse.sh -S /bin/bash -cwd -o $outdir/sge/stdout/parse.out -e $outdir/sge/stderr/parse.err $outdir/sge/scripts/parse.sh`;
