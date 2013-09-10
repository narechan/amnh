#!/usr/bin/perl -w

=head1 NAME

radical_setup.pl

=head1 SYNOPSIS

  radical_setup.pl -- Setup and RADICAL analysis.

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --outdir      Is your output dir
 --config      Is the configuration for your tests
 --step        Is the stepsize across partitions (i.e., 5 to run every 5th partition)
                 (optional)
 --subx        Is set if you want sub-matrices. 

The aln file must be in nexus format.

Dependencies:

Requires the bioperl libraries
Requires Radical.pm

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
use Radical;

my ($help, $matrixfile, $outdir, $configfile, $step, $subx);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    's|step=s'        => \$step,
    'x|subx'          => \$subx,
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
my $cfiobj = Radical->new;
$cfiobj->load_config ($configfile);
$cfiobj->load_aln    ($matrixfile);

# get data we need for sectioned matrices
my $sections = $cfiobj->get_sections;
my $charsets = $cfiobj->get_charsets;
my $mode     = $cfiobj->get_mode;
my $taxa     = $cfiobj->get_taxa;
my $datatype = $cfiobj->get_datatype;
my @charsets = keys (%$charsets);
my $charsetnum = @charsets;

# output the sorted taxa list for later topo analysis
open (T, ">$outdir/summaries/taxa");
foreach my $tax (@$taxa){
    print T "$tax\n";
}
close (T);

# put step to default 1 if not specified
($step = 1) unless ($step);

# generate tree building expts for each sample and
# each partition inclusion
my $sampleMaster = $cfiobj->generate_cmds ($mode,
					   $charsets,
					   $charsetnum,
					   $step,
					   $outdir,
					   $matrixfile);

# generate sub-matrices if required
if ($subx){
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
				 $datatype,
				 "$outdir/sub-matrices");
	}
    }
}
