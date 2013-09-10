#!/usr/bin/perl -w

=head1 NAME

radical_topos.pl

=head1 SYNOPSIS

  radical_topos.pl -- Mine the RADICAL analysis for topologies. Operates
    only on one topology at a time.

Options:

 --help        Show brief help and exit
 --tree        Is your input tree file
 --list        Is a taxa list given in sort order desired for newick serialization

The aln file must be in nexus format.

Dependencies:

Requires the bioperl libraries
Requires Bio::Phylo

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
use Bio::Phylo::IO 'parse';

my ($help, $treefile, $list);
GetOptions(
    'h|help'          => \$help,
    't|tree=s'        => \$treefile,
    'l|list=s'        => \$list,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($treefile, $list){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

my @sort;
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    push (@sort, $line);
}
close (L);

my $proj = parse(
        '-format' => 'nexus',
        '-file'   => $treefile,
        '-as_project' => 1,
                 );

my ($forest)  = @{$proj->get_forests};
my ($tree)    = @{$forest->get_entities};
my $srttree   = $tree->sort_tips(\@sort);
my $srttreest = $srttree->to_newick;

# strip any node labels
$srttreest =~s/\:\d+\.\d+//g;

print "$srttreest\n";
