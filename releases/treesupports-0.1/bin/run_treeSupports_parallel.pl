#!/usr/bin/perl -w

=head1 NAME

run_treeSupports_parallel.pl

=head1 SYNOPSIS

run_treeSupports_parallel.pl

Options:

--prefix       is the SGE run prefix
--matrix       is the data matrix
--tree         is the total evidence tree    
--config       is your paup config
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
--batch        if you want to batch the BS calculations 
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

my ($help, $treefile, $matrixfile, $configfile, $outdir, $kind, $window, $motion, $start, $end, $prefix, $batch);
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
    'p|prefix=s'       => \$prefix,
    'b|batch'        => \$batch,
    ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $matrixfile, $configfile, $outdir, $kind, $prefix){
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
# always done in parallel
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
    
    `qsub -r yes -N $prefix.pbs.$pbsexpt -S /bin/bash -cwd -o $outdir/sge/stdout/$pbsexpt.out -e $outdir/sge/stderr/$pbsexpt.err $outdir/sge/scripts/$pbsexpt.pbs.sh`;
}


# run the run for each BS expt
# BS expts rely on trees from pbs -- so all these must run only
# after PBS runs are complete
opendir (BS, "$outdir/bs/expts");
my @bsexpts = grep (/^.+\..+$/, readdir(BS));
closedir (BS);

# create shell scripts for each BS expt and submit
# as parallel or in batch
if ($batch){
    warn "BS\n";
    open (SHELL, ">$outdir/sge/scripts/bs.sh");
    print SHELL "#!/bin/bash\n";
    print SHELL "date\n";
    print SHELL "source ~/.bashrc\n";

    foreach my $bsexpt (@bsexpts){
        if ( ($window) and ($motion) and ($start) and ($end) ){
            print SHELL "treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x bs -k $kind -o $outdir -f $outdir/bs/expts/$bsexpt -w $window -j $motion -s $start -e $end\n";
        }
        else{
            print SHELL "treeSupports_run.pl -m $matrixfile -t $treefile -c $configfile -x bs -k $kind -o $outdir -f $outdir/bs/expts/$bsexpt\n";
        }
    }
    
    print SHELL "hostname\n";
    print SHELL "date\n";
    close (SHELL);
    
    `qsub -hold_jid $prefix.pbs.* -r yes -N $prefix.bs.sh -S /bin/bash -cwd -o $outdir/sge/stdout/bs.out -e $outdir/sge/stderr/bs.err $outdir/sge/scripts/bs.sh`;
}

else {
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
	
	`qsub -hold_jid $prefix.pbs.* -r yes -N $prefix.bs.$bsexpt -S /bin/bash -cwd -o $outdir/sge/stdout/$bsexpt.out -e $outdir/sge/stderr/$bsexpt.err $outdir/sge/scripts/$bsexpt.bs.sh`;
    }
}

# run the parse
warn "Parse\n";
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

`qsub -hold_jid $prefix.* -r yes -N parse.sh -S /bin/bash -cwd -o $outdir/sge/stdout/parse.out -e $outdir/sge/stderr/parse.err $outdir/sge/scripts/parse.sh`;
