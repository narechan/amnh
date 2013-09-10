#!/usr/bin/perl -w

=head1 NAME

blast_forked.pl

=head1 SYNOPSIS

    This program does a local blast of a fasta file or directory
    of fasta files against a user defined blast database
    and parses the results into a tab-delimited output file.
    
Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta directory or fasta file (script can handle both)
 --blastdb     Is your blastdb
 --config      Is the configuration for your blast
 --outdir      Is your output dir
 --procs       Is the number of forks to run

ncbi blast executables must be in your path
or blastall must be specified in your config

The config must also specify parameters for blast
except the the query and the reference db

BLASTALL=... (optional)
PARAMETERS=...

=head1 DESCRIPTION

Given a dir of fasta files and a  blast db, 
    do and parse the blast.

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

my ($help, $fasta, $blastdb, $config, $outdir, $procs);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'b|blastdb=s'     => \$blastdb,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $blastdb, $config, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/raw`;
`mkdir -p $outdir/parse`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# check if the fasta is a single file or a directory of files
my @queries;
my $signal = 0;
if (-f $fasta){
    push (@queries, $fasta);
    $signal = 1;
}
elsif (-d $fasta){
    opendir (Q, "$fasta");
    @queries = grep (/^.+\..+$/, readdir(Q));
    closedir (Q);
    $signal = 2;
}
else {
    print STDERR "Can't interpret your fasta(s)\n";
    die;
}

# do the forked blast
# note that there is no screen for blast
# against self
my $pm = Parallel::ForkManager->new($procs);
foreach my $query (@queries){
    $pm->start and next;
    
    print STDERR "Blasting $query\n";
    my $name;
    if ($signal == 2){
	$name = blast ($conf, "$fasta/$query", "$blastdb", "$outdir/raw"); 
    }
    elsif ($signal == 1){
	$name = blast ($conf, "$query", "$blastdb", "$outdir/raw");
    }
    
    print STDERR "Parsing $query\n";
    parse ("$outdir/raw", $name, "$outdir/parse");
    
    $pm->finish;
}
$pm->wait_all_children;

#####SUBS#####

sub parse {
    my $indir  = shift;
    my $infile = shift;
    my $outdir = shift;
    my $in = Bio::SearchIO->new(-file   => "$indir/$infile",
				-format => 'blast');
    
    open (OUT, ">$outdir/$infile.parse");
    while( my $result = $in->next_result ) {
        my $hit_count = 0;  
        unless ($result->hits) {
            print OUT join(
			   "\t",
			   $result->query_name,        #1
			   $result->query_length,      #2
			   'No hits found',            #3
			   '0',                        #4
			   '1',                        #5
			   '1',                        #6
			   '1000',                     #7
			   '0',                        #8
			   '0',                        #9
			   '0',                        #10
			   '0',                        #11
			   '0',                        #12
			   '0',                        #13
			   '0',                        #14
			   '0',                        #15
			   '0',                        #16
			   '0',                        #17
			   '0',                        #18
			   '0',                        #19
			   $result->query_description, #20
			   'NULL',                     #21
			   ), "\n";
        }
	
	while( my $hit = $result->next_hit ) {
	    $hit_count++;
	    
	    my $hsp_count = 0;
	    while (my $hsp = $hit->next_hsp ) {
		$hsp_count++;
		print OUT join(
			       "\t",
			       $result->query_name,           #1
			       $result->query_length,         #2
			       $hit->name,                    #3
			       $hit->length(),                #4
			       $hit_count,                    #5
			       $hsp->rank,                    #6
			       $hsp->evalue(),                #7
			       $hsp->score,                   #8
			       $hsp->frac_identical('total'), #9
			       $hsp->start('query'),          #10
			       $hsp->end('query'),            #11
			       $hsp->gaps('query'),           #12
			       $hsp->frac_identical('query'), #13 won't be accurate for blastx
			       $hsp->strand('query'),         #14
			       $hsp->start('hit'),            #15
			       $hsp->end('hit'),              #16
			       $hsp->gaps('hit'),             #17
			       $hsp->frac_identical('hit'),   #18
			       $hsp->strand('hit'),           #19
			       $result->query_description,    #20
			       $hit->description,             #21
			       ), "\n";
		
	    }
	    
	}
	
    }
}

sub blast{
    my $conf   = shift;
    my $fasta  = shift;
    my $db     = shift;
    my $out    = shift;

    my $fasta_name;
    if ($fasta =~/\//g){
	$fasta =~m/.*\/(.*)$/;
	$fasta_name = $1;
    }

    else {
	$fasta_name = $fasta;
    }

    my $db_name;
    if ($db =~/\//g){
	$db =~m/.*\/(.*)$/;
	$db_name = $1;
    }

    else {
	$db_name = $db;
    }

    my $cmd = $conf->{BLASTALL};
    $cmd .= " $conf->{PARAMETERS}";
    $cmd .= " -i $fasta";
    $cmd .= " -d $db";
    $cmd .= " -o $out/$fasta_name.$db_name.out";
    `$cmd`;
    
    return ("$fasta_name.$db_name.out");
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

    # check for blastall in user's path if not in config                                     
    if (! exists($config{"BLASTALL"})){
        my $blastall = `which blastall`;
        chomp $blastall;

        if (defined ($blastall)){
            $config{"BLASTALL"} = $blastall;
        }
        else {
            die "no BLASTALL specified!\n";
	}
    }

    return (\%config);
    close (F);
}
