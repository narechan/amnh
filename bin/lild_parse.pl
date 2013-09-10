#!/usr/bin/perl -w

=head1 NAME

lild_parse.pl

=head1 SYNOPSIS

lild_parse.pl

Options:

--outdir is your output dir
--config is your configfile
--tree is your treefile
--matrix is your alignfile

Requires the bioperl libs. 
Requires Lild.pm

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
use Lild;

my ($help, $treefile, $matrixfile, $configfile, $outdir);
GetOptions(
    'h|help'           => \$help,
    't|tree=s'         => \$treefile,
    'm|matrix=s'       => \$matrixfile,
    'c|config=s'       => \$configfile,
    'o|outdir=s'       => \$outdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($treefile, $matrixfile, $configfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $lildobj = Lild->new;
$lildobj->load_config ($configfile);
$lildobj->load_aln    ($matrixfile);
$lildobj->load_tree   ($treefile);

# cycle through the output files, parse, and
# store the tree lengths
opendir (R, "$outdir/logs");
my @logs = sort (readdir(R));
shift @logs;
shift @logs;
closedir (R);

open (RES, ">$outdir/res/pvals.txt");
foreach my $log (@logs){
    my $data = $lildobj->parse_lild ("$outdir/logs/$log");
    my ($partition, $node, $l) = split (/\./, $log);
    
    open (E, "$outdir/expts/$partition.$node.txt");
    my $coords;
    while (my $line = <E>){
	my @data = split (/\t/, $line);
	$coords = $data[1];
    }
    close (E);

    if (%$data){
    	# calculate length diffence for test part
	my $testdiff = $data->{'TEST'}->{1} -
	    $data->{'TEST'}->{0};
	my $testdiffnorm = $testdiff / $data->{'TEST'}->{0};
	
	# calculate length diffs for rand parts
	my $counter     = 0;
	my $pvalcounter = 0;
	foreach my $part (keys %$data){
	    next if $part eq "TEST";
	    $counter++;
	    
	    my $randdiff = $data->{$part}->{1} -
		$data->{$part}->{0};
	    ($pvalcounter++) if ($randdiff >= $testdiff);
	}
    
	my $pval = sprintf ("%.2f", ($pvalcounter / $counter));
	
	print RES "$partition\t$coords\t$node\t$testdiff\t$testdiffnorm\t";
	print RES "$pvalcounter/$counter\t$pval\n";
    }
    else {
	print RES "$partition\t$coords\t$node\tPAUP errors\n";
    }
    
}
close (RES);
