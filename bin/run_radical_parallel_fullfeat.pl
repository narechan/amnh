#!/usr/bin/perl -w

=head1 NAME

run_radical_parallel.pl

=head1 SYNOPSIS

  run_radical_parallel.pl -- Runs Radical on an SGE cluster.

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --tree        Is your TE treefile (optional)
 --outdir      Is your output dir
 --config      Is the configuration for your tests
 --prefix      Is your sge job name prefix
 --submit      Is your submit strategy 
    (1=atomize every calc; 2=atomize only the samples; 3=hybrid submit; 4=all serial)
    (hybrid submit requires that -x (cutoff) be defined)
    (expts before the cutoff are submitted serially, those after in parallel)
 --cutoff      Is your cutoff (optional)
 --step        Is the stepsize across partitions (i.e., 5 to run every 5th partition)                
                 (optional)
               NOTE THAT IF YOU HAVE A STEP, YOU CAN ONLY USE SUBMIT STRATEGIES 2 AND 4 FOR NOW!
 --virt        Is the amount of memory you require if such a requirement exists (optional)
 --subx        Is set only if you want sub-matrices for each section of the curve
 --doit        Is the list of expts that need doing (works only for serial submit w/o job array)
    (optional)

The aln file must be in nexus format.

Dependencies:
PAUP or RAxML must be installed and in $PATH.                                                                 
Requires Radical.pm                                                                                            
Requires the bioperl libraries                                                                               
Requires Statistics::Descriptive

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

my ($help, $matrixfile, $treefile, $outdir, $configfile, $prefix, $submit, $cutoff, $step, $virt, $subx, $doit);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    't|tree=s'        => \$treefile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    'p|prefix=s'      => \$prefix,
    's|submit=s'      => \$submit,	   
    'x|cutoff=s'      => \$cutoff,
    'y|step=s'        => \$step,
    'v|virt=s'        => \$virt,
    'a|subx'          => \$subx,
    'd|doit=s'            => \$doit,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrixfile, $outdir, $configfile, $prefix, $submit){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

# make sge specific dirs
`mkdir -p $outdir/sge/scripts`;
`mkdir -p $outdir/sge/stdout`;
`mkdir -p $outdir/sge/stderr`;

#####MAIN#####

# run the setup script
# put step to default 1 if not specified                                                                  
=head
($step = 1) unless ($step);
print STDERR "Setup\n";
if ($subx){
   `radical_setup.pl -m $matrixfile -c $configfile -o $outdir -s $step -x`;
}
else{
    `radical_setup.pl -m $matrixfile -c $configfile -o $outdir -s $step`;
}
=cut
# store all the experiments and find reps
opendir (E, "$outdir/cmds");
my @expt = sort (readdir(E));
shift @expt;
shift @expt;
closedir (E);

my %samps;
my %expts;
foreach my $file (@expt){
    my ($samp, $expt, $nex) = split (/\./, $file);
    $samps{$samp} = 1;
    $expts{$expt} = 1;
}

my $samps = keys %samps;
my $expts = keys %expts;
my @expts = sort {$a <=> $b} keys %expts; # using array in case step defined

# create a shell script and submit samples as array (single array of serialized jobs, samples parallelized)
if ($submit == 2){
    open (SHELL, ">$outdir/sge/scripts/submit.sh");
    print SHELL "#!/bin/bash\n";
    print SHELL "date\n";
    print SHELL "source ~/.bashrc\n";

    foreach my $i (@expts){
	print SHELL "radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/\$SGE_TASK_ID.$i.cmds\n";
    }
    print SHELL "hostname\n";
    print SHELL "date\n";
    close (SHELL);
    
    print STDERR "Submitting submit.sh\n";

    if ($virt){
	`qsub -l virtual_free=$virt -r yes -t 1-$samps:1 -N $prefix.cfi.submit.sh -S /bin/bash -cwd -o $outdir/sge/stdout/submit.out -e $outdir/sge/stderr/submit.err $outdir/sge/scripts/submit.sh`;
    }
    else {
	`qsub -r yes -t 1-$samps:1 -N $prefix.cfi.submit.sh -S /bin/bash -cwd -o $outdir/sge/stdout/submit.out -e $outdir/sge/stderr/submit.err $outdir/sge/scripts/submit.sh`;
    }
}

# create a shell script for every samp and submit expts as an array (multiple arrays of parallel jobs)
elsif (($submit == 1) and ($step == 1)){
    for (my $i = 1; $i <= $samps; $i++){
	open (SHELL, ">$outdir/sge/scripts/$i.sh");
	print SHELL "#!/bin/bash\n";
	print SHELL "date\n";
	print SHELL "source ~/.bashrc\n";
	print SHELL "radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$i.\$SGE_TASK_ID.cmds\n";
	print SHELL "hostname\n";
	print SHELL "date\n";
	close (SHELL);
	
	if ($virt){
	    `qsub -l virtual_free=$virt -r yes -t 1-$expts:1 -N $prefix.cfi.$i.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.out -e $outdir/sge/stderr/$i.err $outdir/sge/scripts/$i.sh`;
	}
	else {
	    `qsub -r yes -t 1-$expts:1 -N $prefix.cfi.$i.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.out -e $outdir/sge/stderr/$i.err $outdir/sge/scripts/$i.sh`;
	}
    }
}    

# create a shell script for every sample/number pair and submit sequentially without a job array
# good for if you have a step size
elsif (($submit == 1) and ($step > 1)){
    
    my %exemp;
    open (EX, "$doit");
    while (my $line = <EX>){
	chomp $line;
	$exemp{$line} = 1;
    }
    close (EX);

    for (my $i = 1; $i <= $samps; $i++){
        foreach my $j (@expts){
	    my $string = $i . "." . $j;
	    next unless (exists ($exemp{$string}));

	    open (SHELL, ">$outdir/sge/scripts/$i.$j.sh");
	    print SHELL "#!/bin/bash\n";
	    print SHELL "date\n";
	    print SHELL "source ~/.bashrc\n";
	    print SHELL "time radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$i.$j.cmds\n";
	    print SHELL "hostname\n";
	    print SHELL "date\n";
	    close (SHELL);

	    print STDERR "Submitting $i.$j.sh\n";

	    if ($virt){
		`qsub -l virtual_free=$virt -r yes -N $prefix.cfi.$i.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.$j.out -e $outdir/sge/stderr/$i.$j.err $outdir/sge/scripts/$i.$j.sh`;
	    }
	    else {
		`qsub -r yes -N $prefix.cfi.$i.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.$j.out -e $outdir/sge/stderr/$i.$j.err $outdir/sge/scripts/$i.$j.sh`;
	    }
	}
    }
}

# create a shell script for a single serial submit (no parallelization)
elsif ($submit == 4){
    open (SHELL, ">$outdir/sge/scripts/submit.sh");
    print SHELL "#!/bin/bash\n";
    print SHELL "date\n";
    print SHELL "source ~/.bashrc\n";
    for (my $i = 1; $i <= $samps; $i++){
	foreach my $j (@expts){
	    print SHELL "radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$i.$j.cmds\n";
	}
    }
    print SHELL "hostname\n";
    print SHELL "date\n";
    close (SHELL);

    print STDERR "Submitting submit.sh\n";

    if ($virt){
	`qsub -l virtual_free=$virt -r yes -N $prefix.cfi.submit.sh -S /bin/bash -cwd -o $outdir/sge/stdout/submit.out -e $outdir/sge/stderr/submit.err $outdir/sge/scripts/submit.sh`;
    }
    else {
	`qsub -r yes -N $prefix.cfi.submit.sh -S /bin/bash -cwd -o $outdir/sge/stdout/submit.out -e $outdir/sge/stderr/submit.err $outdir/sge/scripts/submit.sh`;
    }
}


# do a hybrid submit
elsif (($submit == 3) and ($cutoff)){

    # serial portion (parallel on samps)
    open (SHELL, ">$outdir/sge/scripts/submit.sh");
    print SHELL "#!/bin/bash\n";
    print SHELL "date\n";
    print SHELL "source ~/.bashrc\n";

    for (my $i = 1; $i <= $cutoff; $i++){
        print SHELL "radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/\$SGE_TASK_ID.$i.cmds\n";
    }
    print SHELL "hostname\n";
    print SHELL "date\n";
    close (SHELL);

    print STDERR "Submitting submit.sh\n";

    if ($virt){
	`qsub -l virtual_free=$virt -r yes -t 1-$samps:1 -N $prefix.cfi.submit.sh -S /bin/bash -cwd -o $outdir/sge/stdout/submit.out -e $outdir/sge/stderr/submit.err $outdir/sge/scripts/submit.sh`;
    }
    else {
	`qsub -r yes -t 1-$samps:1 -N $prefix.cfi.submit.sh -S /bin/bash -cwd -o $outdir/sge/stdout/submit.out -e $outdir/sge/stderr/submit.err $outdir/sge/scripts/submit.sh`;
    }

    # parallel portion (parallel on samps and expts)
    my $next = $cutoff + 1;
    
    for (my $i = 1; $i <= $samps; $i++){
        open (SHELL2, ">$outdir/sge/scripts/$i.sh");
        print SHELL2 "#!/bin/bash\n";
        print SHELL2 "date\n";
        print SHELL2 "source ~/.bashrc\n";
        print SHELL2 "radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$i.\$SGE_TASK_ID.cmds\n";
        print SHELL2 "hostname\n";
        print SHELL2 "date\n";
        close (SHELL2);

	print STDERR "Submitting $i.sh\n";
	
	if ($virt){
	    `qsub -l virtual_free=$virt -r yes -t $next-$expts:1 -N $prefix.cfi.$i.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.out -e $outdir/sge/stderr/$i.err $outdir/sge/scripts/$i.sh`;
	}
	else {
	    `qsub -r yes -t $next-$expts:1 -N $prefix.cfi.$i.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.out -e $outdir/sge/stderr/$i.err $outdir/sge/scripts/$i.sh`;
	}
    }
}

# or die if submit strategy is unrecognized    
else {
    print STDERR "Unknown submit strategy!\n";
    die;
}

# run the topos scripts
open (T, ">$outdir/sge/scripts/topos.sh");
print T "#!/bin/bash\n";
print T "date\n";
print T "source ~/.bashrc\n";
print T "radical_topos.pl -m $matrixfile -c $configfile -o $outdir\n";
print T "hostname\n";
print T "date\n";
close (T);

print STDERR "Submitting topos.sh\n";
`qsub -hold_jid $prefix.cfi.* -r yes -N topos.sh -S /bin/bash -cwd -o $outdir/sge/stdout/topos.out -e $outdir/sge/stderr/topos.err $outdir/sge/scripts/topos.sh`;

# run the parse when all runs have completed
# need to detect the topology.distribution before continuing
# unless a user input tree is specified                                                              
my $bail = 0;
until ($bail > 0){
    if (-e "$outdir/topology.distribution"){
	sleep 30;
	
	if ($treefile){
	    warn "CFI parse $treefile\n";
	    `radical_parse.pl -m $matrixfile -t $treefile -c $configfile -o $outdir`;
	    $bail++;
	}
	else {
	    my $maxtopos = 0;
	    my %topos;
	    open (TD, "$outdir/topology.distribution");
	    while (my $line = <TD>){
		chomp $line;
		my ($num, $topo) = split (/\t/, $line);
		($maxtopos = $num) if ($num > $maxtopos);
		$topos{$num} = $topo;
	    }
	    close (TD);
	    
	    open (PAR, ">$outdir/sge/scripts/parse.sh");
	    print PAR "#!/bin/bash\n";
	    print PAR "date\n";
	    print PAR "source ~/.bashrc\n";
	    print PAR "radical_parse.pl -m $matrixfile -c $configfile -o $outdir -i \$SGE_TASK_ID -d $outdir/topology.distribution\n";# -p\n";
	    print PAR "hostname\n";
	    print PAR "date\n";
	    close (PAR);
	    
	    warn "Submitting parse.sh\n";
	    `qsub -r yes -t 1-$maxtopos:1 -N parse.sh -S /bin/bash -cwd -o $outdir/sge/stdout/parse.out -e $outdir/sge/stderr/parse.err $outdir/sge/scripts/parse.sh`;
	    $bail++;
	}
    }
    else {
	my $time = time;
	print STDERR "Waiting for topos $time\n";
	sleep 300;
    }
}
