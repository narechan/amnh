#!/usr/bin/perl -w

=head1 NAME

ffp_runner.pl

=head1 SYNOPSIS

  ffp_runner.pl -- wraps around the ffp scripts and creates 
    a distance matrix

Options:

 --help        Show brief help and exit
 --indir       Contains the genomes / texts to be analyzed
 --mode        Is the ffp program used for your datatype
 --length      Is the length of your desired word size
 --outdir      Is the outdir

the ffp related scripts must be in your path.

=head1 DESCRIPTION

Run ffp on files in your indir

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

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $indir, $mode, $length, $outdir);
GetOptions(
    'h|help'          => \$help,
    'i|indir=s'       => \$indir,
    'm|mode=s'        => \$mode,
    'l|length=s'      => \$length,
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;

pod2usage if $help;

for my $option ($indir, $mode, $length, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/vectors`;

#####MAIN#####

opendir (I, "$indir");
my @files = sort (readdir (I));
shift @files;
shift @files;
closedir (I);

foreach my $file (@files){
    warn "FFP $file\n";
    `$mode -l $length $indir/$file > $outdir/vectors/$file.vectors`;
}

opendir (V, "$outdir/vectors");
my @vectors = sort (readdir (V));
shift @vectors;
shift @vectors;
closedir (V);

open (C, ">$outdir/all.vectors");
my @order;
foreach my $vector (@vectors){
    my ($name, $vec) = split (/\./, $vector);
    push (@order, $name);

    open (VF, "$outdir/vectors/$vector");
    while (my $line = <VF>){
	chomp $line;
	print C "$line\n";
    }
    close (VF);
}
close (C);

`ffpcol -t $outdir/all.vectors > $outdir/all.vectors.col`;
`ffprwn $outdir/all.vectors.col > $outdir/all.vectors.row`;
`ffpjsd $outdir/all.vectors.row > $outdir/all.matrix`;

open (NEX, ">$outdir/all.nex");
open (AM, "$outdir/all.matrix");
while (my $line = <AM>){
    chomp $line;
    
    my $file = shift @order;
    print NEX "$file\t$line\n";
}
close (AM);
close (NEX);
