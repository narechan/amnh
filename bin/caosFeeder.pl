#!/usr/bin/perl -w

=head1 NAME

caosFeeder.pl

=head1 SYNOPSIS

  caosFeeder.pl

Options:

 --help        Show brief help and exit
 --treefile    Is your input tree
 --alignfile   Is your input aln
 --outdir      Is your output dir
 --compound    Include if you want to look for
               compound diagostics (warning: mem intense)
 --mapping     Include if you have complex accessions
               and want to run caos through a mapping file

The tree file and aln file must be in nexus format.

CAOS must be in your path.
The treefile and the alignfile MUST have the same accession strings.

=head1 DESCRIPTION

Given an input tree and an input aln nexus file, 
generate CAOS diagnostics and stats.

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
use Bio::AlignIO;

my ($help, $treefile, $alignfile, $outdir, $compound);
GetOptions(
    'h|help'          => \$help,
    't|treefile=s'    => \$treefile,
    'a|alignfile=s'   => \$alignfile,
    'o|outdir=s'      => \$outdir,
    'c|compound'      => \$compound,
    ) or pod2usage;

pod2usage if $help;

#for my $option ($treefile, $alignfile, $outdir){
#    (warn ("Missing a required option\n") and pod2usage)
#	unless ($option);
#}

`mkdir -p $outdir`;

#####MAIN#####

# parse some basic stuff out of the nexus file
warn "Getting some data from your nexus file.\n";
my ($charset, $ntax, $nchar, $datatype, $missing, $gap, $outgroup) = 
    parse_nexus ($alignfile);

# parse the alignment
warn "Parsing your alignment\n";
my ($alignment) = parse_aln ($alignfile, 'nexus');

# parse the tree
#warn "Parsing your tree\n";
#my $tree = parse_tree ($treefile);

# build your caos input file
warn "Building the CAOS input file\n";
open (CAOS, ">$outdir/caos.aln");
#my $taxa = keys (%$alignment);
#print CAOS "$taxa\t$tree\n";
my $counter = 0;
foreach my $taxon (sort keys %$alignment){
    $counter++;
    print CAOS "$counter\t$taxon\t$alignment->{$taxon}\n";
    print "$taxon\n";
}
close (CAOS);    

# run caos with compounds turned on
#warn "Running CAOS\n";
#if ($compound){
#    `CAOS $outdir/caos.aln 2`;
#}
#else {
#    `CAOS $outdir/caos.aln 1`;
#}

#####SUBS#####

sub parse_nexus{
    my $alignfile = shift;
    open (NEX, "$alignfile");
    
    my $ntax;
    my $nchar;
    my $datatype;
    my $missing;
    my $gap;
    my %charset;
    my $outgroup;
    while (my $line = <NEX>){
	chomp $line;
	($nchar    = $1) if ($line =~m/nchar\s*=\s*(\d+)/i);
	($ntax     = $1) if ($line =~m/ntax\s*=\s*(\d+)/i);
	($datatype = $1) if ($line =~m/datatype\s*=\s*(\w+)/i);
	($missing  = $1) if ($line =~m/missing\s*=\s*(.{1})/i);
	($gap      = $1) if ($line =~m/gap\s*=\s*(.{1})/i);

	if ($line =~m/outgroup/i){
	    $line =~s/outgroup//ig;
	    $line =~s/\s+//g;
            $line =~s/\;//g;
	    
	    $outgroup = $line;
	}
	
	
	if ($line =~m/charset/i){
	    $line =~s/charset//ig;
	    $line =~s/\s+//g;
	    $line =~s/\;//g;
	    
	    my ($partition, $coords) = 
		split (/\=/, $line);
	    $charset{$partition} = $coords;
	}
    }
    close (NEX);
    return (\%charset, $ntax, $nchar, $datatype, $missing, $gap, $outgroup);
}

sub parse_aln{
    my $alignfile = shift;
    my $informat  = shift;
    
    my %alignment;
    my $alnin = Bio::AlignIO->new(-file   => "$alignfile",
				  -format => "$informat");

    while (my $aln = $alnin->next_aln()){

	my $counter = 0;
	foreach my $seq ($aln->each_seq){
	    $counter++;

	    # note that the id is stripped of quotes 
	    # and subs spaces for underscores
	    my $id          = $seq->display_id;
	    my $sequence    = $seq->seq;
	    $alignment{$id} = $sequence;
	    
	}
    }
    return (\%alignment);
}

# some of this is taken from orthologID!
sub parse_tree{
    my $treeFile = shift;

    open (FH, "$treeFile");
    my $inTranslate = 0;   # inside translate block                                                          
    my @taxa;
    my $untransTree;
    while (my $line = <FH>) {
	# Lines are not chomp'ed                                                                         
	
	if ($line =~/^\s*Translate/i) {
	    $inTranslate = 1;
	}
	
	elsif ($inTranslate) {
	    my ($key, $val);
	    
	    if ($line =~/^\s*;/){
		$inTranslate = 0; 
		next;
	    }
	    
	    # taken from nexus.pm (TreeIO)
	    if ($line  =~/^\s*(\S+)\s+(.+)$/){
		($key, $val) = ($1, $2);
		$val =~s/[\s,]+$//;
		
	    }
	    
	    $taxa[$key] = $val;
	}
	
	elsif ($line =~/^\s*tree[^\(]+(.*);/i) {
	    ($untransTree = $1) =~ s/\d+/$taxa[$&]/ge;
	}
    }
    return $untransTree;
}
