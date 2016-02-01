#!/usr/bin/perl -w

=head1 NAME

ld_massivePW.pl

=head1 SYNOPSIS

ld_massivePW.pl

Options:

--config       is the config file (treeseach params)
--matrix       is the data matrix
--outdir       is your output dir for the run
--procs       are the number of procs you want to use

the matrix must be in
nexus format.

Requires the bioperl libs.
Requires Lild.pm

=head1 DESCRIPTION

This program sets up one batch ILD calc per gene against all others in charsets.
Using SGE job arrays. So each successive batch run contains less runs. Note that this does not 
run HOMPART, but instead calculates an ILD statistic sans pval.

There is a separate program to run the parse and calc the ILDs

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
use Ild;

my ($help, $configfile, $matrixfile, $outdir, $procs);
GetOptions(
    'h|help'           => \$help,
    'c|config=s'       => \$configfile,   
    'm|matrix=s'       => \$matrixfile,
    'o|outdir=s'       => \$outdir,
    'p|procs=s'       => \$procs,
    ) or pod2usage;
pod2usage if $help;

for my $option ($configfile, $matrixfile, $outdir, $procs){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}


# create dir structure for the results
`mkdir -p $outdir/pw_cmds`;
`mkdir -p $outdir/pw_logs`;
`mkdir -p $outdir/sg_cmds`;
`mkdir -p $outdir/sg_logs`;

#####MAIN#####

# instantiate the object and load required data
my $ildobj = Ild->new;
$ildobj->load_config ($configfile);
$ildobj->load_aln    ($matrixfile);

# get the charsets
my $charsets = $ildobj->get_charsets;

# generate all the single gene commands
print STDERR "Setting up SG analysis\n";
my $sgcounter = 0;
foreach my $charset (sort keys %$charsets){
    $sgcounter++;
    $ildobj->generate_ildstat_singlegene ($sgcounter,
					  $charset,
					  "$outdir/sg_cmds",
					  "$outdir/sg_logs",
					  $matrixfile);
}

# run all the single gene commands
opendir (S, "$outdir/sg_cmds");
my @sgcmds = sort (readdir (S));
shift @sgcmds;
shift @sgcmds;
closedir (S);    

my $sgpm = Parallel::ForkManager->new($procs);
foreach my $sgcmd (@sgcmds){
    $sgpm->start and next;

    print STDERR "Running SG $sgcmd\n";
    $ildobj->run_paup ("$outdir/sg_cmds/$sgcmd");

    $sgpm->finish;
}
$sgpm->wait_all_children;

# for the pairwise all-all case for defined charsets
# does not build pairwise, isolating matices
print STDERR "Setting up PW analysis\n";
my $counter = 0;
foreach my $charset1 (sort keys %$charsets){
    $counter++;
    $ildobj->generate_ildstat_batch_pw ($counter,
					$charset1,
					$charsets,
					"$outdir/pw_cmds",
					"$outdir/pw_logs",
					$matrixfile);
    
    delete $charsets->{$charset1};
    
    # bail if on the last charset (to self only)
    my $numleft = keys %$charsets;
    last if ($numleft == 1);
}

# run all the pairwise commands
opendir (D, "$outdir/pw_cmds");
my @pwcmds = sort (readdir (D));
shift @pwcmds;
shift @pwcmds;
closedir (D);

my $pm = Parallel::ForkManager->new($procs);
foreach my $pwcmd (@pwcmds){
    $pm->start and next;
    
    print STDERR "Running PW $pwcmd\n";
    $ildobj->run_paup ("$outdir/pw_cmds/$pwcmd");
    
    $pm->finish;
}
$pm->wait_all_children;
