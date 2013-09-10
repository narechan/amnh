#!/usr/bin/perl -w

=head1 NAME

chomp_newick.pl

=head1 SYNOPSIS

chomp_newick.pl

Options:

--treefile is your newick treefile

=head1 DESCRIPTION

This program removes all return characters
from a newick tree

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

my ($help, $treefile);
GetOptions(
    'h|help'          => \$help,
    't|tree=s'        => \$treefile,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($treefile){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

my @taxa;
open (T, "$treefile");
while (my $line = <T>){
    chomp $line;
    $line =~s/\s//g;
    
    if ($line =~m/\'(.*)\'/){
	push (@taxa, $1);
	
    }

    print "$line";
}
print "\n";
