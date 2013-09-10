#!/usr/bin/perl -w

=head1 NAME

missingTaxaMatrixer.pl

=head1 SYNOPSIS

missingTaxaMatrixer.pl

Options:

--matrix is your alignfile
--cutoff is the missing data threshold

The matrix must be in nexus format.

Requires the bioperl libs. 

=head1 DESCRIPTION

This program generates matrices given missing data cutoffs.
Assumes that missing characters are denoted by '?'.

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

my ($help, $matrixfile, $cutoff);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'c|cutoff=s'      => \$cutoff,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $cutoff){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# instantiate the object and load stuff we need
my $supportobj = TreeSupports->new;
$supportobj->load_aln    ($matrixfile);

# get the charsets
my $charsets = $supportobj->get_charsets;
my $chars    = $supportobj->get_nchar;

# store alignment information
my $partitions = {};
my $lengths    = {};
my $alnin = Bio::AlignIO->new(-file   => "$matrixfile",
			      -format => "nexus");

# only one aln there
my $aln = $alnin->next_aln();
my %taxa;
foreach my $seq ($aln->each_seq){
    my $id        = $seq->display_id;
    $taxa{$id} = 1;

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

# cycle through the characters and charsets and count missing
my %freq;
foreach my $cset (sort keys %$partitions){
    my $len    = 0;
    $freq{$cset} = 0;

    foreach my $id (sort keys %{ $partitions->{$cset} }){
	my $seq = $partitions->{$cset}->{$id};
	$len = length ($seq);

	my $matchstring = '\\' . '?' . "{$len}"; # a bit hokey!
	($freq{$cset}++) unless ($seq =~m/$matchstring/);
    }
}

# find the total number of characters given the cutoff
# and build the new partition boundaries
my $newparts = {};
my $start = 1;
my $nchars = 0;
my $taxa = keys %taxa;
foreach my $part (sort keys %freq){
    my $ratio = $freq{$part} / $taxa;
    next if ($ratio < $cutoff);
    
    my $end = $start + ($lengths->{$part} - 1);
    $newparts->{$part} = "$start-$end";
    $nchars += $lengths->{$part};
    $start = $end + 1;
}    

# cycle through the charsets and build new matrices
print "#NEXUS\n";
print "BEGIN DATA;\n";
print "DIMENSIONS NTAX=$supportobj->{'alignment'}->{'ntax'} NCHAR=$nchars;\n";
print "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print "MATRIX\n\n";

foreach my $part (sort keys %freq){
    my $ratio = $freq{$part} / $taxa;
    next if ($ratio < $cutoff);

    print "[Partition $part chars $newparts->{$part}]\n";
    foreach my $tax (sort keys %{ $partitions->{$part} }){
	print "$tax\t$partitions->{$part}->{$tax}\n";
    }
    print "\n";
}
print ";\n";
print "END;\n\n";

print "BEGIN SETS;\n";
foreach my $part (sort keys %$newparts){
    print "CHARSET $part = $newparts->{$part};\n";
}
print "END;\n";

    
