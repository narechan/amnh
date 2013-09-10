#!/usr/bin/perl -w

=head1 NAME

metasim_forked.pl

=head1 SYNOPSIS

  metasim_forked.pl -- this program will generate sequence reads
    using metasim.

Options:

 --help        Show brief help and exit
 --fasta       Is your input genome
 --reads       Is a file that contains the number of reads desired in each
 --config      Is the configuration for your meatsim
 --outdir      Is your output dir
 --procs       Is the number of forks to run

metasim must be in your path

The config must specify parameters for metasim

METASIM=...

=head1 DESCRIPTION

Simulate reads

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

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Bio::SearchIO;

my ($help, $fasta, $reads, $config, $outdir, $procs);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'r|reads=s'       => \$reads,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $reads, $config, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/reads`;
`mkdir -p $outdir/stderr`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# do the forked metasim runs for every read count
my $pm = Parallel::ForkManager->new($procs);
open (R, "$reads");
while (my $line = <R>){
    chomp $line;
    warn "$line $fasta\n";
    
    $pm->start and next;
    metasim ($line, $fasta, $conf, $outdir);
    $pm->finish;
    
}
$pm->wait_all_children;
close (R);

#####SUBS#####

sub metasim{
    my $reads  = shift;
    my $fasta  = shift;
    my $conf   = shift;
    my $out    = shift;

    # make tmp outdir
    `mkdir $out/reads/$reads-tmp`;

    # metasim
    my $metasim = "MetaSim cmd";
    ($metasim .= " $conf->{'METASIM'}") if ($conf->{'METASIM'});
    $metasim .= " -r $reads";
    $metasim .= " -d $out/reads/$reads-tmp";
    $metasim .= " $fasta";
    `$metasim &> $out/stderr/stderr.$reads`;
    `mv $out/reads/$reads-tmp/* $out/reads/$reads.fasta`;
    `rm -rf $out/reads/$reads-tmp/`;
    
}

sub parse_config{
    my $file = shift;

    my %config;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }

    return (\%config);
    close (F);
}
