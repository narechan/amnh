#!/usr/bin/perl -w

=head1 NAME

insidR.pl

=head1 SYNOPSIS

This script parses the PRC summary file from insid.pl, executes the mm.R
R script which models the PRC saturation behavior using a hyperbolic 
equation (the Michaelis-Menten equation), and calculates the query sequence
basepairs and the amount of the reference aligned (or covered) at the 
desired slope or desired query sequence coverage. Tabulated results printed to stdout.

Dependencies:                                                                                            

Requires R.

Options:

--infile is the summary infile from insid.pl
--slope is the slope you want to solve for (optional)
--coverage is the coverage (in total bp count) you want to solve for (optional)
--rscript is the R script you want to run (mm.R: included in the package)
--outdir is the dir where R modelling data will be saved

Note that either slope or coverage must be defined.

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

my ($help, $infile, $slope, $rscript, $outdir, $seqcoverage);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'm|slope=s'      => \$slope,
    'r|rscript=s'     => \$rscript,
    'o|outdir=s'       => \$outdir,
    'c|coverage=s'     => \$seqcoverage,
    ) or pod2usage;
pod2usage if $help;

for my $option ($infile, $rscript, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####    

# get summary file name
my $infile_name;
if ($infile =~/\//g){
    $infile =~m/.*\/(.*)$/;
    $infile_name = $1;
}

else {
    $infile_name = $infile;
}

# run the R script for your summary file
print STDERR "Running R script on $infile\n";
`R --slave --args $infile < $rscript > $outdir/$infile_name.Rout`;

# parse the summary output
my $vmax;
my $km;
open (F, "$outdir/$infile_name.Rout");
while (my $line = <F>){
    chomp $line;

    my @line;
    if ($line =~m/^Vm\s*/){
	@line = split (/\s+/ , $line);
	$vmax = $line[1];
    }
    
    if ($line =~m/^K\s*/){
	@line = split (/\s+/, $line);
	$km = $line[1];
    }
}
close (F);

if ($slope){
    
    # calculate the coverage and RefCovered at the desired slope
    my $coverage = sqrt ((($vmax * $km) / $slope)) - $km;
    my $refcovered = ($vmax * $coverage) / ($km + $coverage);
    print "$infile_name\t$km\t$vmax\t$slope\t$coverage\t$refcovered\n";

}

elsif ($seqcoverage){
    
    # calculate the ref covered and the slope at the desired sequence coverage
    my $refcovered = ($vmax * $seqcoverage) / ($km + $seqcoverage);
    my $slope      = ($vmax * $km) / (($km + $seqcoverage) ** 2);
    print "$infile_name\t$km\t$vmax\t$slope\t$seqcoverage\t$refcovered\n";
    
}

else {
    print STDERR "Need a slope or a sequence coverage.\n";
    die;
}
