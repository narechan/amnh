#!/usr/bin/perl -w

=head1 NAME

insidR.pl

=head1 SYNOPSIS

The script parses the PRC summary file from insid.pl, executes the mm.R
R script which models the PRC saturation behavior using a hyperbolic 
equation (the Michaelis-Menten equation), and calculates the query sequence
basepairs and the amount of the reference aligned (or covered) at the 
desired slope. Tabulated results printed to stdout.

Dependencies:                                                                                            

Requires R.

Options:

--infile/-i is the summary infile from insid.pl
--slope/-m is the slope you want to solve for
--rscript/-r is the R script you want to run (mm.R: included in the package)
--outdir/-o is the dir where R modelling data will be saved

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

my ($help, $infile, $slope, $rscript, $outdir);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'm|slope=s'      => \$slope,
    'r|rscript=s'     => \$rscript,
    'o|outdir=s'       => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($infile, $slope, $rscript, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####    

# find the WGA refcovered from the summary file
# if available
open (W, "$infile");
my $rc;
while (my $line = <W>){
    chomp $line;
    my @line = split (/\t/, $line);
    $rc = $line[3];
}
close (W);

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

# calculate the coverage and RefCovered at the desired slope
my $coverage = sqrt ((($vmax * $km) / $slope)) - $km;
my $refcovered = ($vmax * $coverage) / ($km + $coverage);

# calculate the coverage and slope at the refCovered in the WGA
unless ($rc eq "UNKNOWN"){
    my $coverage1 = ($km * $rc) / ($vmax - $rc);
    my $slope1 = ($vmax * $km) / (($km + $coverage1) ** 2);
    print "$infile_name\t$km\t$vmax\t$slope\t$coverage\t$refcovered\t$slope1\t$coverage1\t$rc\n";
}
else {
    print "$infile_name\t$km\t$vmax\t$slope\t$coverage\t$refcovered\tUNKOWN\tUNKNOWN\tUNKNOWN\n";
}
