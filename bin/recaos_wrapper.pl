#!/usr/bin/perl -w

=head1 NAME

recaos_wrapper.pl

=head1 SYNOPSIS

  recaos_wrapper.pl

Options:

 --help        Show brief help and exit
 --treefile    Is your input tree
 --alignfile   Is your input aln
 --outdir      Is your output dir

The tree file must be in nexus
The aln file must be in fasta

reCAOS must be in your path.
The treefile and the alignfile MUST have the same accession strings.

=head1 DESCRIPTION

Given an input tree and an input aln nexus file, 
generate reCAOS diagnostics.

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

#####TODO:
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::TreeIO;

my ($help, $treefile, $alignfile, $outdir);
GetOptions(
    'h|help'          => \$help,
    't|treefile=s'    => \$treefile,
    'a|alignfile=s'   => \$alignfile,
    'o|outdir=s'      => \$outdir,
    ) or pod2usage;

pod2usage if $help;

for my $option ($treefile, $alignfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# convert the nexus file to a newick string
my $in  = Bio::TreeIO->new(-file  => $treefile, -format => "nexus");
my $out = Bio::TreeIO->new(-format => "newick", -file => ">$outdir/tmpAPN.tre");

my $t = $in->next_tree;
$out->write_tree($t);

# read the alignment, the newick tree,
# and generate reCAOS input file
open (RC, ">$outdir/recaos.in");
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$alignfile");
my $counter = 0;
while (my $sequence_obj = $seqin->next_seq()){
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();

    print RC "$id $seq\n";
}

open (NW, "$outdir/tmpAPN.tre");
while (my $line = <NW>){
    chomp $line;
    chop $line;
    print RC "$line\n"; #expecting just one line from the newick file
}
close (NW);
close (RC);

# run reCAOS
`reCAOS $outdir/recaos.in > $outdir/diag.chars`;
`rm $outdir/tmpAPN.tre`;
