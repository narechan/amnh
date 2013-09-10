#!/usr/bin/perl -w

=head1 NAME

insidResampGrid.pl

=head1 SYNOPSIS

The script wraps around insidP.pl, insid.pl, insidR.pl, and insidD.pl
to either bootstrap or jackknife an insid reads tree. It iterates over the 
entire insid pipeline for every replicate and generates a majority 
rules bootstrap or jackknife tree for the entire process. It operates
in the context of an SGE cluster, parallelizing the resampling of reads
and the insid calculations

Note that for resamp, query coverages must be provided to as total base
counts.

Dependencies:                                                                                            

All insid components and mummer must be in your path
Requires the bioperl libraries.                                                                               
Requires the runner module Mummer.pm.
Requires R.

Options:

--data is a directory containing read fasta files for
    all species
--mode is the sampling mode (bootstrap = 1; jackknife = 2)
--config is the mummer configuration file
--refbases is a file containing the base counts                                                            
    ref subset
--qrybases is a file containing the base counts
    for the qry subset
--outdir is the output dir for all data
--reps is the number of replicates you want
--clean cleans out the mummer aln files as generated (optional)

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2010 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $datadir, $outdir, $refbasefile, $qrybasefile, $mode, $configfile, $reps, $clean);
GetOptions(
    'h|help'          => \$help,
    'd|datadir=s'     =>\$datadir,
    'o|outdir=s'      => \$outdir,
    'a|refbases=s'       => \$refbasefile,
    'b|qrybases=s'    => \$qrybasefile,	   
    'c|config=s'      => \$configfile,
    'm|mode=s'        => \$mode,
    'r|reps=s'        => \$reps,
    'x|clean'         => \$clean,
    ) or pod2usage;
pod2usage if $help;

for my $option ($datadir, $outdir, $refbasefile, $qrybasefile, $mode, $configfile, $reps){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir/reps`;

#####MAIN#####    

# get all orgs
opendir (D, "$datadir");
my @readfiles = sort (readdir (D));
shift @readfiles;
shift @readfiles;
closedir (D);

# get the basecounts
my @refbasecounts;
open (B, "$refbasefile");
while (my $basecount = <B>){
    chomp $basecount;
    push (@refbasecounts, $basecount);
}
close (B);

my @qrybasecounts;
open (Q, "$qrybasefile");
while (my $basecount = <Q>){
    chomp $basecount;
    push (@qrybasecounts, $basecount);
}
close (Q);

for (my $i = 1; $i <= $reps; $i++){
    `mkdir -p $outdir/reps/$i/results`;
    `mkdir -p $outdir/reps/$i/models`;
    `mkdir -p $outdir/reps/$i/matrices`;
    `mkdir -p $outdir/reps/$i/sge/scripts`;
    `mkdir -p $outdir/reps/$i/sge/stderr`;
    `mkdir -p $outdir/reps/$i/sge/stdout`;

    # resample to create query and ref reads sets for this rep (serial submit)
    foreach my $readfile (@readfiles){
	open (RESAMP, ">$outdir/reps/$i/sge/scripts/insidP.$readfile.sh");
	print RESAMP "#!/bin/bash\n";
	print RESAMP "source ~/.bashrc\n";
	print RESAMP "insidP.pl -i $datadir/$readfile -o $outdir/reps/$i/qry_reads/$readfile/ -b $qrybasefile -m $mode\n";
	print RESAMP "insidP.pl -i $datadir/$readfile -o $outdir/reps/$i/ref_reads/$readfile/ -b $refbasefile -m $mode\n";
	close (RESAMP);
	
	print STDERR "Submitting INSIDP for $i $readfile\n";
	`qsub -r yes -N insidP.$i.$readfile -S /bin/bash -cwd -o $outdir/reps/$i/sge/stdout/insidP.$readfile.out -e $outdir/reps/$i/sge/stderr/insidP.$readfile.err $outdir/reps/$i/sge/scripts/insidP.$readfile.sh`;
    }
    
    # build the insid analysis for every set of query reads against every 
    # reference basecount across all pairwise comparisons in this rep.
    # a bash script is created for each analysis.
    my $insidcounter = 0;
    foreach my $basecountR (@refbasecounts){
	foreach my $reffile (@readfiles){
	    foreach my $qryfile (@readfiles){
		$insidcounter++;

		open (INSID, ">$outdir/reps/$i/sge/scripts/insid.$insidcounter.sh");
		print INSID "#!/bin/bash\n";
		print INSID "source ~/.bashrc\n";
		
		print INSID "insid.pl -p 1 -c $configfile -s $outdir/reps/$i/qry_reads/$qryfile -r $outdir/reps/$i/ref_reads/$reffile/$basecountR.fa -o $outdir/reps/$i/results/$basecountR/$qryfile-$reffile\n";
		
		
		print INSID "mv $outdir/reps/$i/results/$basecountR/$qryfile-$reffile/summary $outdir/reps/$i/results/$basecountR/$qryfile-$reffile/$qryfile-$reffile.summary\n";
		close (INSID);
	    }
	}
    }
    
    # submit the insid analysis as an array job
    open (SBMT, ">$outdir/reps/$i/sge/scripts/insid.sh");
    print SBMT "#!/bin/bash\n";
    print SBMT "source ~/.bashrc\n";
    print SBMT "chmod 744 $outdir/reps/$i/sge/scripts/insid.\$SGE_TASK_ID.sh\n";
    print SBMT "./$outdir/reps/$i/sge/scripts/insid.\$SGE_TASK_ID.sh\n";
    close (SBMT);

    print STDERR "Submitting INSID for $i\n";
    `qsub -hold_jid insidP.* -t 1-$insidcounter:1 -r yes -N insid.$i -S /bin/bash -cwd -o $outdir/reps/$i/sge/stdout/insid.out -e $outdir/reps/$i/sge/stderr/insid.out $outdir/reps/$i/sge/scripts/insid.sh`;
 
    # build the matrix analysis for each reference basecount (serial submit)
    foreach my $basecountR (@refbasecounts){
	
	opendir (D, "$outdir/reps/$i/results/$basecountR");
        my @dirs = sort (readdir (D));
        shift @dirs;
        shift @dirs;
        closedir (D);

        my $refcovdata = {};
        my $actualbasecountQs = {};
        foreach my $dir (@dirs){

            open (F, "$outdir/reps/$i/results/$basecountR/$dir/$dir.summary");
            while (my $line = <F>){
                chomp $line;

                my ($reads, $qbases, $rbases, $rcovered, $rfrac) =
                    split (/\t/, $line);
                $refcovdata->{$dir}->{$qbases} = $line;
                $actualbasecountQs->{$qbases} = 1;
            }
            close (F);
        }
	
	# map actual basecounts to input basecounts                     
        my $realQs = {};
        foreach my $basecountQ (@qrybasecounts){
            my $mindiff = 1e99;
            my $realbasecountQ;

            foreach my $actualbasecountQ (sort keys %$actualbasecountQs){
                my $diff = abs($basecountQ - $actualbasecountQ);

                if ($diff < $mindiff){
                    $mindiff = $diff;
                    $realbasecountQ = $actualbasecountQ;
                }
                else {
                    next;
                }
            }
            $realQs->{$basecountQ} = $realbasecountQ;
        }


	foreach my $basecountQ (@qrybasecounts){
	    
	    open (M, ">$outdir/reps/$i/models/$basecountQ-$basecountR.models");
	    foreach my $pwcomp (keys %$refcovdata){
                my $mindiff = 1e99;
                my $realbasecountQ;

                # do the actual basecount mapping                       
                foreach my $actualbasecountQ (keys %{$refcovdata->{$pwcomp}}){
                    my $diff = abs($basecountQ - $actualbasecountQ);

                    if ($diff < $mindiff){
                        $mindiff = $diff;
                        $realbasecountQ = $actualbasecountQ;
                    }
                    else {
                        next;
                    }
                }

                print M "$pwcomp\t$refcovdata->{$pwcomp}->{$realbasecountQ}\n";
            }
            close (M);

	    # build the insidD command
	    open (MATX, ">$outdir/reps/$i/sge/scripts/insidD.$basecountQ-$basecountR.sh");
	    print MATX "#!/bin/bash\n";
	    print MATX "source ~/.bashrc\n";
	    print MATX "insidD.pl -i $outdir/reps/$i/models/$basecountQ-$basecountR.models -o $outdir/reps/$i/matrices/$basecountQ-$basecountR\n";
	    close (MATX);

	    print STDERR "Submitting INSIDD for $i $basecountQ-$basecountR\n";
	    `qsub -hold_jid insid.* -r yes -N insidD.$i.$basecountQ-$basecountR -S /bin/bash -cwd -o $outdir/reps/$i/sge/stdout/insidD.out -e $outdir/reps/$i/sge/stderr/insidD.err $outdir/reps/$i/sge/scripts/insidD.$basecountQ-$basecountR.sh`;
	}
    }

    # stall until insidDs complete
    while (1){
	my $counter = 0;
	my $treecounter = 0;
	foreach my $basecountR (@refbasecounts){
	    foreach my $basecountQ (@qrybasecounts){
		$counter++;
		
		if (-e "$outdir/reps/$i/matrices/$basecountQ-$basecountR/tree"){
		    $treecounter++;
		}
		else {
		    next;
		}
	    }
	}
	
	if ($treecounter == $counter){
	    (`rm -rf $outdir/reps/$i/qry_reads/*`) if ($clean);
	    (`rm -rf $outdir/reps/$i/ref_reads/*`) if ($clean);
	    (`rm -rf $outdir/reps/$i/results/*`) if ($clean);
	    last;
	}
	else {
	    my $time = time;
	    print STDERR "Waiting for rep $i to complete $time\n";
	    sleep 60;
	}
    }
}
