#!/usr/bin/perl -w

=head1 NAME

ddh_wrapper.pl

=head1 SYNOPSIS

ddh_wrapper.pl

Options:

--procs is the number of || processes for ddh.pl
--indir is the dir that contains your input query genomes
--refdir is the dir that contains your read references
--outdir is your output dir

Requires the bioperl libs. 
Requires Mummer.pm, Metasim.pm, Blast.pm, ddh.pl
Requires blast, formatdb, mummer and metasim
    to be in your path

NOTE THAT THERE ARE HARD CODED BITS HERE:
    THE NUMBER OF FORKED DDH.PL'S
    THE NEXT GEN SEQ STRATEGIES

=head1 DESCRIPTION

ddh.pl runs for every ref/query combo of genomes in the indir
once per seq strat

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
use Parallel::ForkManager;

my ($help, $procs, $indir, $outdir, $refdir);
GetOptions(
    'h|help'          => \$help,
    'p|procs=s'       => \$procs,
    'i|indir=s'       => \$indir,
    'o|outdir=s'      => \$outdir,
    'r|refdir=s'      => \$refdir,
    ) or pod2usage;
pod2usage if $help;

for my $option ($procs, $indir, $outdir, $refdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# get the query genomes
opendir (Q, "$indir");
my @genomes = grep (/^.+\..+$/, readdir(Q));
closedir (Q);

# get the reference read sets
opendir (R, "$refdir");
my @refs = sort (readdir (R));
shift @refs;
shift @refs;
closedir (R);

# build a ddh.pl command list
my @commands;
foreach my $qgenome (sort @genomes){
    foreach my $rgenome (sort @refs){
	
	# need to do this so metasim doesn't get confused with punctuation in the dir
	my ($rgenome_name, $rcov) = split (/\./, $rgenome);
	my ($qgenome_name, $qfas) = split (/\./, $qgenome);
	
	if ($qgenome_name eq $rgenome_name){
	    next;
	}
	else {
	    my $command = "ddh_mummer.pl -p $procs -c config.ddh.mummer -q $indir/$qgenome -r $refdir/$rgenome -o $outdir/$qgenome_name-v-$rgenome_name-$rcov >& $outdir/stderr.$qgenome_name-v-$rgenome_name-$rcov"; ##SOME HARD##
	    push (@commands, $command);
	}
    }
}

# run the accumulated commands
my $pm = Parallel::ForkManager->new(1); ##HARD##
foreach my $cmd (@commands){
    $pm->start and next;
#    print STDERR "$cmd\n";
    `$cmd`;
    $pm->finish;
}
$pm->wait_all_children;
