#!/usr/bin/perl -w

=head1 NAME

jk_parse.pl

=head1 SYNOPSIS

jk_parse.pl

Options:

--outdir is your output dir
--config is your configfile
--matrix is your alignfile

The matrix must be in nexus format.
Requires the bioperl libs. 
Requires Jk.pm

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
use Jk;

my ($help, $matrixfile, $configfile, $outdir);
GetOptions(
    'h|help'           => \$help,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $configfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $jkobj = Jk->new;
$jkobj->load_config ($configfile);
$jkobj->load_aln ($matrixfile);

my $tax = $jkobj->get_tax;
my $taxwant = $tax - 3;

# cycle through the output  files, parse, and
# store the ild
opendir (R, "$outdir/logs");
my @logs = sort (readdir(R));
shift @logs;
shift @logs;
closedir (R);

open (RES, ">$outdir/summary");
foreach my $log (@logs){
    my ($pctdel, $loggie) = split (/\./, $log);
    my $freqs = $jkobj->parse_jk ("$outdir/logs/$log");
    
    my $freqsum = 0;
    my $freqcount = 0;
    my $counter = 0;
    foreach my $freq (@$freqs){
	$counter++;
	$freqcount++;
	$freqsum += $freq;
	last if ($counter == $taxwant);
    }
    my $freqav = $freqsum / $freqcount;
    
    my $freqstring = join ",", @$freqs;
    print RES "$pctdel\t$freqstring\t$freqav\n";
}
close (RES);
