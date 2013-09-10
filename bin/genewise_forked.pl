#!/usr/bin/perl -w

=head1 NAME

genewise_forked.pl

=head1 SYNOPSIS

  genewise_forked.pl -- 
              

Options:

 --help        Show brief help and exit
 --fastadb     Is your query fasta directory
    (individual sequences where you want to call genes)
 --aafile      Is your proteins database
    (all proteins homologs/orthologs you want to search for)
 --config      Is the configuration for your genewise
 --outdir      Is your output dir
 --procs       Is the number of forks to run

genewise must be in your path.
make sure genewiseconfigdir is defined in your env.

The config must specify parameters for genewise

PARAMETERS=...

=head1 DESCRIPTION

Given a dir of dna seq and a dir of proteins, 
    find genes using genewise.

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
use Bio::Tools::Genewise;

my ($help, $fastadb, $aafile, $config, $outdir, $procs);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fastadb,
    'o|outdir=s'      => \$outdir,
    'a|aafile=s'      => \$aafile,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fastadb, $outdir, $aafile, $config, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# get the dna seqs
opendir (F, "$fastadb");
my @fastadb = sort (readdir (F));
shift @fastadb;
shift @fastadb;
closedir (F);


# do the forked gene search
=head
my $pm = Parallel::ForkManager->new($procs);
foreach my $fasta (@fastadb){
    $pm->start and next;
    
    print STDERR "Genewise $fasta\n";
    run_genewise ($conf, $aafile, $fastadb, $fasta, $outdir);
	
    $pm->finish;
}
$pm->wait_all_children;
=cut

# do the parse
foreach my $fasta (@fastadb){
    parse_genewise ($outdir, $fasta);
}

#####SUBS#####

sub parse_genewise{
    my $out   = shift;
    my $fasta = shift;
    
    my $gw = Bio::Tools::Genewise->new(-file=>"$out/genewise.out");
    while (my $gene = $gw->next_prediction){
	my @transcripts = $gene->transcripts;
	foreach my $t(@transcripts){
	    my @exons =  $t->exons;
	    foreach my $e(@exons){
		print $e->start." ".$e->end."\n";
	    }
	}
    }
}
    
sub run_genewise{
    my $conf  = shift;
    my $aafile = shift;
    my $fastadb = shift;
    my $fasta = shift;
    my $out   = shift;

    my $cmd = "genewisedb";
    $cmd .= " $conf->{PARAMETERS}";
    $cmd .= " -prodb $aafile";
    $cmd .= " -dnas $fastadb/$fasta";
    $cmd .= " > $out/$fasta.out";
    `$cmd`;
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
    close (F);

    return (\%config);
}
