#!/usr/bin/perl -w

=head1 NAME

insidContig_parallel.pl

=head1 SYNOPSIS

The script wraps around insid.pl and insidD.pl to create a 
matrix and distance tree for pairwise comparisons of contigs

It doesn't really do anything parallel but can be easily submitted
in a parallel manner.

Dependencies:                                                                                            

All insid components and mummer must be in your path
Requires the bioperl libraries.                                                                               
Requires the runner module Mummer.pm.

Options:

--ref is the reference sequence
--data is a directory containing contig fasta files for
    all species
--config is the mummer configuration file
--outdir is the output dir for all data

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
use Parallel::ForkManager;

my ($help, $datadir, $outdir, $configfile, $ref);
GetOptions(
    'h|help'          => \$help,
    'd|datadir=s'     =>\$datadir,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    'r|ref=s'       => \$ref,
    ) or pod2usage;
pod2usage if $help;

for my $option ($datadir, $outdir, $configfile, $ref){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir/results`;
`mkdir -p $outdir/matrices`;

#####MAIN#####    

# get all orgs
opendir (D, "$datadir");
my @readfiles = sort (readdir (D));
shift @readfiles;
shift @readfiles;
closedir (D);

# do insid run for all against reference
foreach my $qryfile (@readfiles){
	
    # insid
    print STDERR "INSID for $reffile $qryfile\n"; 
    `insid.pl -p 1 -t -c $configfile -s $datadir/$qryfile -r $ref -o $outdir/results/$qryfile-$reffile`;
    
    open (O, ">$outdir/results/$qryfile-$reffile/summary.mod"); 
    open (F, "$outdir/results/$qryfile-$reffile/summary");
    my $line = <F>;
    chomp $line;
    my @line = split (/\t/, $line);
    unshift (@line, "$qryfile-$reffile");
    my $newline = join "\t", @line;
    print O "$newline\n";
    close (F);
    close (O);
    
    `rm $outdir/results/$qryfile-$reffile/summary`;
    
    # move the result to the parent directory to deal 
    # with the too many links error (too many subdirs)
    `mv $outdir/results/$qryfile-$reffile/summary.mod $outdir/results/$qryfile-$reffile.summary`;
    `rm -rf $outdir/results/$qryfile-$reffile`;
    
}
