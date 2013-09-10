#!/usr/bin/perl -w

=head1 NAME

submatrixer.pl

=head1 SYNOPSIS

  submatrixer.pl -- generate a submatrix with only the desired partitions

Options:

 --help        Show brief help and exit
 --matrix      Is your input matrix (with partitions defined)
 --golist      Is the list of GO terms and the associated linking gene
 --oglist      Is the list of OGs and their associcated linking genes and go terms   
 --outdir      Is your output directory

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

my ($help, $matrixfile, $outdir, $golistfile, $oglistfile);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'g|golist=s'      => \$golistfile,
    'p|oglist=s'      => \$oglistfile,
    'o|outdir=s'      => \$outdir,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrixfile){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# instantiate the object and load required data                                                       
my $cfiobj = Radical->new;
$cfiobj->load_aln    ($matrixfile);

# get data we need for sectioned matrices
my $charsets = $cfiobj->get_charsets;
my @charsets = keys (%$charsets);
my $charsetnum = @charsets;

# associate linker genes with their OGs
my $linkers = {};
open (L, "$oglistfile");
while (my $line = <L>){
    chomp $line;
    my ($og, $at, $gos) = split (/\t/, $line);
    $linkers->{$at} = $og;
}
close (L);

my @list;
my $counter = 0;
open (G, "$golistfile");
while (my $line = <G>){
    $counter++;
    chomp $line;
    
    my @list;
    my ($goterm, $ats) = split (/\t/, $line);
    $goterm =~s/\:/\_/;
    
    my @ats = split (/\;/, $ats);
    foreach my $at (@ats){
	print STDERR "$goterm\t$at\t$linkers->{$at}\n";
	push (@list, $linkers->{$at});
    }

    my $line = join " ", @list;
    
    $cfiobj->sub_matrix ($goterm,
			 $counter,
			 $line,
			 $matrixfile,
			 $charsets,
			 "$outdir");
    
}
close (G);
