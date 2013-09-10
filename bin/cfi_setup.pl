#!/usr/bin/perl -w

=head1 NAME

cfi_setup.pl

=head1 SYNOPSIS

  run_cfi.pl

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --outdir      Is your output dir
 --config      Is the configuration for your tests
 --step        Is the stepsize across partitions (i.e., 5 to run every 5th partition)
                 (optional)
 --subx        Is set if you want sub-matrices. 
 --cmddir     Set only if you want to define sub-matrices for an existing analysis 
                 (using its cmddir) post-hoc.

The aln file must be in nexus format.

PAUP must be in your path
Requires Cfi.pm

The config must specify all parameters for the
search (treecommand), the reps (samples), 
the max trees searched (maxtrees).

Maxtrees can be left undefined in which case it will be
set to 100 (default) and the command, "SET increase=auto" 
will be used.
Examp:

TREECOMMAND=hsearch swap=nni addseq=closest
MAXTREES=200
SAMPLES=10

=head1 DESCRIPTION

Calculate CFI given successive additions of all partitions,
one at a time

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2009 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::AlignIO;
use Cfi;

my ($help, $matrixfile, $outdir, $configfile, $step, $subx, $cmddir);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    's|step=s'        => \$step,
    'x|subx'          => \$subx,
    'd|cmddir=s'      => \$cmddir,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrixfile, $outdir, $configfile){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/cmds`;
`mkdir -p $outdir/logs`;
`mkdir -p $outdir/trees`;
`mkdir -p $outdir/summaries`;
`mkdir -p $outdir/output`;
`mkdir -p $outdir/consensus`;
(`mkdir -p $outdir/sub-matrices`) if ($subx);

#####MAIN#####

# instantiate the object and load required data                                                       
my $cfiobj = Cfi->new;
$cfiobj->load_config ($configfile);
$cfiobj->load_aln    ($matrixfile);

# get data we need for sectioned matrices
my $sections = $cfiobj->get_sections;
my $charsets = $cfiobj->get_charsets;
my $mode     = $cfiobj->get_mode;
my @charsets = keys (%$charsets);
my $charsetnum = @charsets;

# put step to default 1 if not specified
($step = 1) unless ($step);

# generate tree building expts for each sample and
# each partition inclusion or harvest the nodes used 
# in a prior analysis to regenerate sub-matrices
my $sampleMaster = {};
if ($cmddir){
    opendir (CMDS, "$cmddir");
    my @cmdfiles = sort (readdir (CMDS));
    shift @cmdfiles;
    shift @cmdfiles;
    closedir (CMDS);
    
    foreach my $cmdfile (@cmdfiles){
	my ($sample, $number, $stuff) = split (/\./, $cmdfile);

	open (F, "$cmddir/$cmdfile");
	while (my $line = <F>){
	    chomp $line;
	    $sampleMaster->{$sample}->{$number} = $line;
	}
	close (F);
    }
}
else{
    $sampleMaster = $cfiobj->generate_cmds ($mode,
					    $charsets,
					    $charsetnum,
					    $step,
					    $outdir,
					    $matrixfile);
}

# generate sub-matrices if required
if ($subx){

    # get the ends of each section                                                                            
    my %ends;
    if ($sections =~m/\,/){
	$sections =~s/\s//g;
	my @sections = split (/\,/, $sections);
	foreach my $sect (@sections){
	    $ends{$sect} = 1;
	}
    }
    else {
	my $sectparts = round ($charsetnum / $sections);
	my $start = 1;
	my $end = 0;
	until ($start > $charsetnum){
	    $end = $start + $sectparts - 1;
	    if ($end > $charsetnum){
		$end = $charsetnum;
	    }
	    
	    $ends{$end} = 1;
	    $start = $end + 1;
	}
    }
    
    foreach my $sample (sort {$a <=> $b} keys %$sampleMaster){
        foreach my $number (sort {$a <=> $b} keys %{$sampleMaster->{$sample}}){
            next unless (exists ($ends{$number}));

	    $cfiobj->sub_matrix ($sample,
				 $number,
				 $sampleMaster->{$sample}->{$number},
				 $matrixfile,
				 $charsets,
				 "$outdir/sub-matrices");
	}
    }
}
