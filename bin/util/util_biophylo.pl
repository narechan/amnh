#!/usr/bin/perl
use strict;
use warnings;
use Bio::Phylo::IO 'parse';

my $file = shift @ARGV;
my $proj = parse(
        '-format' => 'nexus',
        '-file'   => $file,
        '-as_project' => 1,
		 );

my ($forest) = @{ $proj->get_forests };
my ($tree) = @{ $forest->get_entities };
my $ladder = $tree->ladderize;
print $file, "\t", $ladder->to_newick, "\n";
