#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Phylo::IO 'parse';

my ($help, $treedir);
GetOptions(
    'h|help'             => \$help,
    't|treedir=s'        => \$treedir,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($treedir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# read in tree files
opendir (D, "$treedir");
my @treefiles = sort (readdir(D));
shift @treefiles;
shift @treefiles;
closedir (D);

# convert tree files to newick and consolidate
# strings into a single file
open (T, ">/tmp/trees.newick");
foreach my $tree (@treefiles){
    my $proj = parse(
		     '-format' => 'nexus',
		     '-file'   => "$treedir/$tree",
		     '-as_project' => 1,
                 );

    my ($forest)  = @{$proj->get_forests};
    my ($tree)    = @{$forest->get_entities};
    my $treest    = $tree->to_newick;

    print T "$treest\n";
}
close (T);

# do consensus of the newick topologies
my $forest = parse(
		   -format => 'newick',
		   -file   => "/tmp/trees.newick",
		   );

my $consensus = $forest->make_consensus( -branches => 'frequency' );
my $cnstreest = $consensus->to_newick;
print "$cnstreest\n";

