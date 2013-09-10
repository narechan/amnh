#!/usr/bin/perl -w

=head1 NAME

radical_check.pl

=head1 SYNOPSIS

  radical_check.pl -- check if the radical library satisfies the statistical
    criterion for abort at the concat point where instantiated.

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --outdir      Is your output dir
 --config      Is the configuration for your tests
 --tree        Is the topology to check the replicates against
 --point       Is the concatenation point you want to asses

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
use Statistics::Descriptive;

my ($help, $matrixfile, $outdir, $configfile, $treefile, $point);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
    't|treefile=s'    => \$treefile,
    'p|point=s'       => \$point,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrixfile, $outdir, $configfile, $treefile, $point){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# instantiate the object and load required data                                                          
my $cfiobj = Radical->new;
$cfiobj->load_config ($configfile);
$cfiobj->load_aln    ($matrixfile);
my $root = $cfiobj->get_root;
my $taxa = $cfiobj->get_taxa;
my $abortFrac = $cfiobj->get_abortFrac;
my $abortInt  = $cfiobj->get_abortInt;

# get all the trees                                                                                       
opendir (TREES, "$outdir/trees");
my @trees = sort (readdir (TREES));
shift @trees;
shift @trees;
closedir (TREES);

# parse the TE tree and get max cfi
my $topoTE = $cfiobj->parse_tree ($treefile, $root, "nexus");
my $maxcfi = @$taxa - 3;

# get the CFI stats across all replicates for 
# the concat point of interest
open (S, ">>$outdir/check.sensor");
my $samps = 0;
my $maxsamps = 0;
my @cfis;
foreach my $tree (@trees){
    my ($rep, $concatpt, $nex) = split (/\./, $tree);
    next unless ($concatpt eq $point);
    $samps++;
    
    # get the topo of the rep tree
    my $topoREP = $cfiobj->parse_tree ("$outdir/trees/$tree", $root, "nexus");

    # calculate the CFI
    my $cfidata;
    my ($t, $topocomp) = $cfiobj->compare_all_nodes ($topoTE, $topoREP);
    foreach my $top (sort keys %$t){
	$cfidata += $t->{$top};
    }
    $cfidata--;
    $cfidata--;
    push (@cfis, $cfidata);
    
    ($maxsamps++) if ($cfidata == $maxcfi);
}    

my $ratio = $maxsamps / $samps;
my $cfis = join "\t", @cfis;
print S "CURRENTCFIs:$cfis\n";
print S "MAXCFIs:$point\t$ratio\n";
close (S);

# run the abort sensor
my $points;
open (C, "$outdir/check.sensor");
while (my $line = <C>){
    chomp $line;
    if ($line =~m/MAXCFIs/){
	$line =~s/MAXCFIs\://;
	my ($pt, $ratio) = split (/\t/, $line);
	$points->{$pt} = $ratio;
    }
    else {
	next;
    }
}
close (C);
 
my $status;
my $statuscnt = 0;
foreach my $pt (sort {$a <=> $b} keys %$points){
    if ($points->{$pt} >= $abortFrac){
	$status = "stable";
	$statuscnt++;
    }
    else{
	$status = "unstable";
    }
}

open (X, ">>$outdir/check.sensor");
if (($status eq "stable") and ($statuscnt >= $abortInt)){
    print X "STABLE\n";
}
else {
    print X "UNSTABLE\n";
}
close (X);
