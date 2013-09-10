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
 --submit      Is your submit strategy (1=atomize every calc; 2=atomize only the replicates)
 --step        Is the stepsize across partitions (i.e., 5 to run every 5th partition)                
                 (optional)
 --virt        Is the amount of memory you require if such a requirement exists (optional)
 --subx        Is set only if you want sub-matrices for each section of the curve
 --checking    Is set if you want to abort concatenation when the statistical criterion 
    has been reached for the replicates. If you want to use the abort procedure, a treefile 
    must be specified.

The aln file must be in nexus format.

Dependencies:
PAUP or RAxML must be installed and in $PATH.                                                                 
Requires Radical.pm                                                                                            
Requires the bioperl libraries                                                                               
Requires Statistics::Descriptive
Requires Bio::Phylo

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

my ($help, $matrixfile, $treefile, $outdir, $configfile, $prefix, $step, $virt, $subx, $checking, $submit);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    't|tree=s'        => \$treefile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    'p|prefix=s'      => \$prefix,
    'v|virt=s'        => \$virt,
    'a|subx'          => \$subx,
    'z|checking'      => \$checking, 	   
    's|submit=s'      => \$submit,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrixfile, $outdir, $configfile, $prefix){
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
($step = 1) unless ($step);
print STDERR "Setup\n";
if ($subx){
   `radical_setup.pl -m $matrixfile -c $configfile -o $outdir -s $step -x`;
}
else{
    `radical_setup.pl -m $matrixfile -c $configfile -o $outdir -s $step`;
}

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

# submit samples
my $concats = {};
if ($checking){
    open (C, ">$outdir/sge/scripts/checker.sh");
    print C "#!/bin/bash\n";
    print C "date\n";
    print C "source ~/.bashrc\n";
    print C "radical_check.pl -m $matrixfile -t $treefile -c $configfile -o $outdir > $outdir/temp.stats\n";
    print C "hostname\n";
    print C "date\n";
    close (C);

    foreach my $j (@expts){
	$concats->{$j} = 1;
	
	# checking after every parallel submit of a concatenation point (samples parallelized)
	if ($submit == 1){
	    foreach my $i (sort {$a <=> $b} keys %samps){
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
		    `qsub -l virtual_free=$virt -r yes -N $prefix.$i.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.$j.out -e $outdir/sge/stderr/$i.$j.err $outdir/sge/scripts/$i.$j.sh`;
		}
		else {
		    `qsub -r yes -N $prefix.$i.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.$j.out -e $outdir/sge/stderr/$i.$j.err $outdir/sge/scripts/$i.$j.sh`;
		}
	    }
	}

	# checking after every serial submit of a concatenation point
	elsif ($submit == 2){
	    open (SHELL, ">$outdir/sge/scripts/$j.sh");
	    print SHELL "#!/bin/bash\n";
	    print SHELL "date\n";
	    print SHELL "source ~/.bashrc\n";

	    foreach my $i (sort {$a <=> $b} keys %samps){
		print SHELL "time radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$i.$j.cmds\n";
	    }
	    print SHELL "hostname\n";
	    print SHELL "date\n";
	    close (SHELL);
	    
	    print STDERR "Submitting $j.sh\n";
	    
	    if ($virt){
		`qsub -l virtual_free=$virt -r yes -N $prefix.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$j.out -e $outdir/sge/stderr/$j.err $outdir/sge/scripts/$j.sh`;
	    }
	    else {
		`qsub -r yes -N $prefix.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$j.out -e $outdir/sge/stderr/$j.err $outdir/sge/scripts/$j.sh`;
	    }
	}
	else {
	    print STDERR "Submit strategy not recognized\n";
	    die;
	}

	# submit the checker
	print STDERR "Submitting checker.sh\n";
	`qsub -hold_jid $prefix.* -r yes -N checker.sh -S /bin/bash -cwd -o $outdir/sge/stdout/checker.out -e $outdir/sge/stderr/checker.err $outdir/sge/scripts/checker.sh`;
	
	# hold until temp.stats appears for this iteration
	my $bail = 0;
	until ($bail > 0){
	    if (-e "$outdir/temp.stats"){
		sleep 30;
		$bail++;
	    }
	    else {
		sleep 30;
		next;
	    }
	}
	
	# check if the tree lib has stabilized
	open (S, "$outdir/temp.stats");
        my $state = <S>;
        chomp $state;
        warn "$state\n";

        if ($state =~m/^STABLE/){
            `rm $outdir/temp.stats`;
            last;
        }
        else {
            `rm $outdir/temp.stats`;
            next;
        }
	
    }
}
else {

    # parallelize all replicates and concatenation points in a complete radical without checking
    if ($submit == 1){
	foreach my $j (@expts){
	    $concats->{$j} = 1;
	    
	    foreach my $i (sort {$a <=> $b} keys %samps){
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
		    `qsub -l virtual_free=$virt -r yes -N $prefix.$i.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.$j.out -e $outdir/sge/stderr/$i.$j.err $outdir/sge/scripts/$i.$j.sh`;
		}
		else {
		    `qsub -r yes -N $prefix.$i.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.$j.out -e $outdir/sge/stderr/$i.$j.err $outdir/sge/scripts/$i.$j.sh`;
		}
	    }
	}
    }

    # parallelize only the replicates in a complete radical without checking
    elsif ($submit == 2){
	foreach my $i (sort {$a <=> $b} keys %samps){
	    open (SHELL, ">$outdir/sge/scripts/$i.sh");
	    print SHELL "#!/bin/bash\n";
	    print SHELL "date\n";
	    print SHELL "source ~/.bashrc\n";
	    
	    foreach my $j (@expts){
		$concats->{$j} = 1;
		print SHELL "time radical_run.pl -m $matrixfile -c $configfile -o $outdir -e $outdir/cmds/$i.$j.cmds\n";
	    }
	    print SHELL "hostname\n";
	    print SHELL "date\n";
	    close (SHELL);
	    
	    print STDERR "Submitting $i.sh\n";
	    
	    if ($virt){
		`qsub -l virtual_free=$virt -r yes -N $prefix.$i.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.out -e $outdir/sge/stderr/$i.err $outdir/sge/scripts/$i.sh`;
	    }
	    else {
		`qsub -r yes -N $prefix.$i.sh -S /bin/bash -cwd -o $outdir/sge/stdout/$i.out -e $outdir/sge/stderr/$i.err $outdir/sge/scripts/$i.sh`;
		
	    }
	}
    }
    else {
	print STDERR "Submit strategy not recognized\n";
	die;
    }
}

# run the topos scripts on each concat pt serially
my $treecounter = 0;
foreach my $j (sort {$a <=> $b} keys %$concats){
    open (TOPOS, ">$outdir/sge/scripts/topos.$j.sh");
    print TOPOS "#!/bin/bash\n";
    print TOPOS "date\n";
    print TOPOS "source ~/.bashrc\n";
    
    foreach my $i (sort {$a <=> $b} keys %samps){
	$treecounter++;
	print TOPOS "time radical_topos.pl -l $outdir/summaries/taxa -t $outdir/trees/$i.$j.tre > $outdir/$i.$j.temp.tre\n";
    }
    
    print TOPOS "hostname\n";
    print TOPOS "date\n";
    close (TOPOS);
    
    print STDERR "Submitting topos for $j\n";
    
    `qsub -hold_jid $prefix.* -r yes -N topos.$j.sh -S /bin/bash -cwd -o $outdir/sge/stdout/topos.$j.out -e $outdir/sge/stderr/topos.$j.err $outdir/sge/scripts/topos.$j.sh`;
}

# check to see whether all topos have been written
# when topos set, concatenate, and begin the post process
my $baily = 0;
until ($baily > 0){
    my $tracker = 0;
    foreach my $j (sort {$a <=> $b} keys %$concats){
	foreach my $i (sort {$a <=> $b} keys %samps){
	    if (-e "$outdir/$i.$j.temp.tre"){
		$tracker++;
	    }
	    else {
		next;
	    }
	}
    }
    
    if ($tracker == $treecounter){
	sleep 30;
	
	# create the topology distribution
	my $topologies = {};
	opendir (DIR, "$outdir");
	my @files = sort (readdir (DIR));
	shift @files;
	shift @files;
	closedir (DIR);
	
	foreach my $file (@files){
	    next unless ($file =~m/temp\.tre/);
	    
	    open (TF, "$outdir/$file");
	    while (my $line = <TF>){
		chomp $line;
		$line =~s/\;//;
		$topologies->{$line}++;
	    }
	    close (TF);
	    `rm $outdir/$file`;
	}

	open (TS, ">$outdir/topology.distribution");
	my $tcounter = 0;
	foreach my $topo (keys %$topologies){
	    $tcounter++;
	    print TS "$tcounter\t$topo\t$topologies->{$topo}\n";
	}
	close (TS);

	# run the parse
	if ($treefile){
	    warn "CFI post-process $treefile\n";
	    `radical_pp.pl -m $matrixfile -t $treefile -c $configfile -o $outdir`;
	    $baily++;
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
	    
	    open (PAR, ">$outdir/sge/scripts/pp.sh");
	    print PAR "#!/bin/bash\n";
	    print PAR "date\n";
	    print PAR "source ~/.bashrc\n";
	    print PAR "radical_pp.pl -m $matrixfile -c $configfile -o $outdir -i \$SGE_TASK_ID -d $outdir/topology.distribution\n";# -p\n";
	    print PAR "hostname\n";
	    print PAR "date\n";
	    close (PAR);
	    
	    warn "Submitting pp.sh\n";
	    `qsub -r yes -t 1-$maxtopos:1 -N parse.sh -S /bin/bash -cwd -o $outdir/sge/stdout/pp.out -e $outdir/sge/stderr/pp.err $outdir/sge/scripts/pp.sh`;
	    $baily++;
	}
    }
    else {
	my $time = time;
	print STDERR "Waiting for topos $time\n";
	sleep 60;
    }
}
