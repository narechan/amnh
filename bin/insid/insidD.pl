#!/usr/bin/perl -w

=head1 NAME

insidD.pl

=head1 SYNOPSIS

The script creates a distance matrix using the COGWES values generated 
by insidR.pl. It also generates the UPGMA tree of this distance matrix 
and outputs a tabulated distance table.

Dependencies:                                                                                            

None.

Options:

--infile is the summary infile from insidR.pl for reads
    or insid.pl for contigs, concatenated over all pairwise comparisons
--outdir is the directory for output files
--contigs is set if you're analyzing contig data directly from insid.pl

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

my ($help, $infile, $outdir, $contigs);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'o|outdir=s'      => \$outdir,
    't|contigs'       => \$contigs,
    ) or pod2usage;
pod2usage if $help;

for my $option ($infile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####    

# parse and store the data
my $data = {};
open (F, "$infile");
while (my $line = <F>){
    chomp $line;
    
    my @data = split (/\t/, $line);
    $data[0]=~s/\.summary//;
    my ($reference, $query) = split (/\-/, $data[0]);
    
    if ($contigs){
	$data->{$reference}->{$query} = $data[3];
    }
    else{
	$data->{$reference}->{$query} = $data[5];
    }
}
close (F);

# get the number of taxa
my @taxa = keys %$data;
my $taxa = @taxa;

# harvest self comparisons for reference levels
my $reflevels = {};
foreach my $r (sort (keys %$data)){
    foreach my $q (sort (keys %{$data->{$r}})){
	($reflevels->{$r} = $data->{$r}->{$q}) if 
	    ($r eq $q);
   }
}

# build the distance matrix
open (M, ">$outdir/matrix");
open (T, ">$outdir/table");
print M "#NEXUS\n";
print M "Begin distances;\n";
print M "Dimensions ntax=$taxa;\n";
print M "format triangle=both;\n";
print M "matrix\n";

my $counter = 0;
foreach my $ref (sort (keys %$data)){
    $counter++;
    print M "$ref\t";

    my @cogdists;
    foreach my $qry (sort (keys %{$data->{$ref}})){

	my $qrycog = $data->{$qry}->{$ref};
	my $refcog = $data->{$ref}->{$qry};
	
	my $avgcog = ($qrycog + $refcog) / 2;
	my $normfactor = ($reflevels->{$ref} + $reflevels->{$qry}) / 2;
	print STDERR "$qry\t$ref\t$qrycog\t$refcog\t$avgcog\t$normfactor\n";
	my $normavgcog = $avgcog / $normfactor;
	my $cogdist = sprintf ("%.15f", (1 - $normavgcog));
	
	push (@cogdists, $cogdist);
	print T "$qry\t$ref\t$cogdist\n";
    }

    my $cogdiststring = join "\t", @cogdists;

    if ($counter == $taxa){
	print M "$cogdiststring;\n";
    }
    else {
	print M "$cogdiststring\n";
    }
}

print M "end;\n";
print M "Begin paup;\n";
print M "UPGMA;\n";
print M "savetree file=$outdir/tree;\n";
print M "end;\n";

`paup -n $outdir/matrix`;
