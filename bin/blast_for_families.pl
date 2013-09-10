#!/usr/bin/perl -w

=head1 NAME

blast_for_families.pl

=head1 SYNOPSIS

  blast_for_families.pl 
              

Options:

 --help        Show brief help and exit
 --fasta       Is your fasta file of candidate genes
 --blastdbs    Is your the directory containing all your blastdbs
 --config      Is the configuration for your blast
 --outdir      Is your output dir
 --eval        Is your eval cutoff for reporting results
 --genomes     Is the directory of genomes/assemblies

ncbi blast executables must be in your path
or blastall must be specified in your config

The config must also specify parameters for blast
except the the query and the reference db

BLASTALL=... (optional)
PARAMETERS=...

=head1 DESCRIPTION

Given a fasta file and a collection of genomes or assemblies,
blast the fasta against the assemblies and generate putative 
orthologous groups by extracting sequence from source fasta files.
Parses a TILING PATH from the hsps.

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
use Bio::SearchIO;
use Bio::Index::Fasta;
use Bio::Location::Simple;

my ($help, $fasta, $blastdbs, $config, $outdir, $evalcutoff, $genomes);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'b|blastdbs=s'    => \$blastdbs,
    'c|config=s'      => \$config,
    'e|evalcutoff=s'  => \$evalcutoff,
    'g|genomes=s'     => \$genomes,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $blastdbs, $config, $evalcutoff, $genomes){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/blast/raw`;
`mkdir -p $outdir/blast/parse`;
`mkdir -p $outdir/families`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# get all the blastdbs
my %bdbs;
opendir (Q, "$blastdbs");
my @blastdbs = grep (/^.+\..+$/, readdir(Q));
foreach my $d (@blastdbs){
    my @name = split (/\./, $d);
    pop @name;
    
    my $b = join ".", @name;
    $bdbs{$b} = 1;
}
closedir (Q);

# do the blast against all dbs
my $families = {};
foreach my $db (sort (keys %bdbs)){
    
    print STDERR "Blasting $fasta $db\n";
    my $name = blast ($conf, "$fasta", "$blastdbs/$db", "$outdir/blast/raw"); 
        
    print STDERR "Parsing $fasta $db\n";
    my $parse = parse ("$outdir/blast/raw", $name, "$outdir/blast/parse", $evalcutoff, $db);
    
    # transpose into $families
    foreach my $db (keys %$parse){
	foreach my $gene (keys %{$parse->{$db}}){
	    $families->{$gene}->{$db} = $parse->{$db}->{$gene};
	}
    }
}

# extract the sequence from the genomes to 
# create family fasta files
foreach my $gene (keys %$families){
    open (GF, ">$outdir/families/$gene.fa");
    
    foreach my $db (keys %{$families->{$gene}}){
	next if (($families->{$gene}->{$db} eq "NO HITS") or
                 ($families->{$gene}->{$db} eq "LOW HITS"));
	
	my $index = Bio::Index::Fasta->new(-filename => "$genomes/$db" . ".idx", 
					   -write_flag => 1);
	$index->make_index("$genomes/$db");
	
	my $hit = $families->{$gene}->{$db}[0];
	my $start = $families->{$gene}->{$db}[1];
	my $end = $families->{$gene}->{$db}[2];
	my $ori = $families->{$gene}->{$db}[3];
	
	my $location = Bio::Location::Simple->new(-start  => $start,
                                                  -end    => $end,
                                                  -strand => $ori);
	my $sequence = $index->fetch($hit);
        my $subseq   = $sequence->subseq($location);
	print GF ">$db\n$subseq\n";
#	`rm $genomes/$db
    }
}

#####SUBS#####

sub parse {
    my $indir  = shift;
    my $infile = shift;
    my $outdir = shift;
    my $evalcutoff = shift;
    my $db = shift;

    my $in = Bio::SearchIO->new(-file   => "$indir/$infile",
				-format => 'blast');

    my $parse = {};
    open (F, ">$outdir/$infile.parse");
    while (my $result = $in->next_result){
	my $queryname = $result->query_name;
	my $querylen  = $result->query_length;
	
	# bail if there are no hits
	unless ($result->hits){
	    print F "$queryname\tNo hits found\n";
	    $parse->{$db}->{$queryname} = 'NO HITS';
	    next;
	}

	# get only the best hit
	my @hits = $result->hits;
	my $besthit = shift @hits;
	my $hitname = $besthit->name;
	
	# cycle through the hsps and 
	# confirm that at least one exceeds the eval threshold.
	# if so, compile and report tiled coordinates
	my $minstart = 1000000000000000000000000000000;
	my $maxend   = 0;
	my $mineval  = 100000000000000000;
	my $ori;
	while (my $hsp = $besthit->next_hsp){
	    my $evalue = $hsp->evalue;
	    my $hspstart = $hsp->start('hit');
	    my $hspend   = $hsp->end('hit');
	    my $hspori   = $hsp->strand('hit');
	    
	    ($minstart = $hspstart) if ($hspstart < $minstart);
	    ($maxend   = $hspend) if ($hspend > $maxend);
	    ($mineval  = $evalue) if ($evalue < $mineval);
	    $ori = $hspori;
	}
	
	my $hsplen = $maxend - $minstart + 1;
	my $lenreq = $hsplen / $querylen;
	
	# report if the evalue cutoff satisfied at least once
	# and within length tolerances
	if (($mineval < $evalcutoff) and ($lenreq > 0.80) and ($lenreq < 1.20)){
	    print F "$queryname\t$hitname\t$minstart\t$maxend\t$ori\n";
	    
	    my @data;
	    push (@data, $hitname, $minstart, $maxend, $ori);
	    $parse->{$db}->{$queryname} = \@data;
	}
	else {
	    print F "$queryname\t$hitname\tDoes not satisfy hit criteria\n";

	    $parse->{$db}->{$queryname} = 'LOW HITS'; 
	}
    }
    close (F);
    
    return ($parse);
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
