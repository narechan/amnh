#!/usr/bin/perl -w

=head1 NAME

rmChars.pl

=head1 SYNOPSIS

rmChars.pl

Options:

--matrix is your alignfile
--missing is your missing char

The matrix must be in nexus format.

Requires the bioperl libs. 

=head1 DESCRIPTION

This program finds all charsets that have full
sequence representation.

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
use TreeSupports;
use Bio::AlignIO;
use Bio::SeqIO;

my ($help, $matrixfile, $missing);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'n|missing=s'     => \$missing,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $missing){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $supportobj = TreeSupports->new;
$supportobj->load_aln    ($matrixfile);

# get the charsets
my $charsets = $supportobj->get_charsets;

# store alignment information
my $partitions = {};
my $lengths    = {};
my $alnin = Bio::AlignIO->new(-file   => "$matrixfile",
			      -format => "nexus");

# only one aln there
my $aln = $alnin->next_aln();
foreach my $seq ($aln->each_seq){
    my $id        = $seq->display_id;
    
    foreach my $charset (sort keys %$charsets){
#	print STDERR "Storing\t$id\t$charset\n";

	my $coords = $charsets->{$charset};
	my ($start, $end) = split (/\-/, $coords);
	
	my $partition = $seq->subseq($start, $end);
	my $partlen   = length ($partition);
	
	$partitions->{$charset}->{$id} = $partition;
	$lengths->{$charset} = $partlen;
    }
}

# cycle through the charsets and see which 
# have data for all their taxa
my @completes;
foreach my $cset (sort keys %$partitions){

    my $signal = 0;
    foreach my $id (sort keys %{ $partitions->{$cset} }){
	my $seq = $partitions->{$cset}->{$id};
	my $len = length ($seq);
	
	my $matchstring = '\\' . $missing . "{$len}"; # a bit hokey!

	($signal = 1) if ($seq =~m/$matchstring/);
    }
 
    (print "$cset\n") unless ($signal == 1);
}
