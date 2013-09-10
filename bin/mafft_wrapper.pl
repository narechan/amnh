#!/usr/bin/perl -w

=head1 NAME

mafft_wrapper.pl

=head1 SYNOPSIS

  mafft_wrapper.pl

Options:

 --help        Show brief help and exit
 --indir       Contains the files to be analyzed
 --config      Is the configuration for your run (MAFFT params)
 --outdir      Is your output dir
 --procs       Is the number of procs you want to use

mafft must be in your path

=head1 DESCRIPTION

Usage examp:

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
use Bio::SeqIO;
use Bio::AlignIO;
use Parallel::ForkManager;

my ($help, $indir, $config, $outdir, $procs);
GetOptions(
    'h|help'          => \$help,
    'i|indir=s'       => \$indir,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    ) or pod2usage;

pod2usage if $help;

for my $option ($indir, $outdir, $config, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/alns`;
`mkdir -p $outdir/stderr`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# do the mafft alignments
opendir (A, "$indir");
my @fastas = sort readdir (A);
shift @fastas;
shift @fastas;
closedir (A);
my $pm = Parallel::ForkManager->new($procs);
foreach my $fasta (sort @fastas){
    $pm->start and next;
    
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$indir/$fasta");
    my $seqcount = 0;
    while (my $sequence_obj = $seqin->next_seq()){
	$seqcount++;
    }
    if ($seqcount > 1){
	print STDERR "Mafft $fasta\n";
	mafft_run ($fasta, $conf, $outdir, $indir);
    }
    else {
	print STDERR "Mafft skipped $fasta\n";
    }
    $pm->finish;
}
$pm->wait_all_children;


#####SUBS#####

sub mafft_run{
    my $fas  = shift;
    my $conf = shift;
    my $out  = shift;
    my $in   = shift;

    # mafft
    my $mafft = "mafft";
    ($mafft .= " $conf->{'MAFFT'}") if ($conf->{'MAFFT'});
    $mafft .= " $in/$fas";

    `$mafft 1>$out/alns/$fas.aln 2>$out/stderr/$fas.stderr`;
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

