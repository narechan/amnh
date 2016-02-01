#!/usr/bin/perl -w

=head1 NAME

missingTaxaAnalyzer.pl

=head1 SYNOPSIS

missingTaxaAnalyzer.pl

Options:

--matrix is your alignfile
--outdir is your outdir

The matrix must be in nexus format.

Requires the bioperl libs. 

=head1 DESCRIPTION

This program quantifies missing data across the entire matrix and
per charset. Assumes that missing characters are denoted by '?'.

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

my ($help, $matrixfile, $outdir);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'o|outdir=s'      => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

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
foreach my $seq ($aln->each_seq){
    my $id        = $seq->display_id;
    
    foreach my $charset (sort keys %$charsets){
	print STDERR "Storing\t$id\t$charset\n";

	my $coords = $charsets->{$charset};
	my ($start, $end) = split (/\-/, $coords);
	
	my $partition = $seq->subseq($start, $end);
	my $partlen   = length ($partition);
	
	$partitions->{$charset}->{$id} = $partition;
	$lengths->{$charset} = $partlen;
    }
}

# cycle through the characters and charsets and count missing
open (PARTS, ">$outdir/partition.stats");
my %matrix;
my %freq;
my $fulls;
foreach my $cset (sort keys %$partitions){
    my $coords = $charsets->{$cset};
    my ($start, $end) = split (/\-/, $coords);

    my $len    = 0;
    $freq{$cset} = 0;
    
    my $fulltaxa = 0;
    my @fulltaxaids;

    open (F, ">$outdir/$cset.fa");
    foreach my $id (sort keys %{ $partitions->{$cset} }){
	my $seq = $partitions->{$cset}->{$id};
	$len = length ($seq);

	print F ">$id\n$seq\n";

	# full taxa per partition
	my $matchstring = '\\' . '?' . "{$len}"; # a bit hokey!                   
	($freq{$cset}++) unless ($seq =~m/$matchstring/);

	# per strain counts
	($fulltaxa++) unless ($seq =~m/$matchstring/);
#	($fulltaxaids->{$id} = 1) unless ($seq =~m/$matchstring/);
	(push (@fulltaxaids, $id)) unless ($seq =~m/$matchstring/);

	# partition level counts
	my $missingCount = ($seq =~ tr/?//);
	print PARTS "$cset\t$id\t$missingCount\n";

	# matrix level counts
	my $counter = $start;
	for (split //, $seq){
	    ($matrix{$counter}++) if ($_ eq '?');
	    $counter++;
	}
	
    }
    close (F);

    # per strain data collation
    if ($fulltaxa <= 10){
	foreach my $fullt (@fulltaxaids){
	    $fulls->{$fullt}->{$cset} = $fulltaxa;
	}
    }
    else {
	next;
    }
}
close (PARTS);

open (PTAXA, ">$outdir/pertaxa.stats");
foreach my $taxa (keys %$fulls){
    my $pcounter = 0;
    my @string;
    foreach my $gene (keys %{$fulls->{$taxa}}){
	$pcounter++;
	my $string = $gene . ":" . $fulls->{$taxa}->{$gene};
	push (@string, $string);
    }
    my $jstring = join ",", @string;
    print PTAXA "$taxa\t$pcounter\t$jstring\n";
}
close (PTAXA);

open (ALL, ">$outdir/matrix.stats");
#foreach my $position (sort {$a <=> $b} keys %matrix){
for (my $i = 1; $i <= $chars; $i++){
    if (exists ($matrix{$i})){
	print ALL "$i\t$matrix{$i}\n";
    }
    else {
	print ALL "$i\t0\n";
    }
}
close (ALL);

open (F, ">$outdir/freqs.stats");
foreach my $freq (sort keys %freq){
    print F "$freq\t$freq{$freq}\n";
}
close (F);

