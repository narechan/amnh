#!/usr/bin/perl -w

=head1 NAME

blast_with_ref_cov_forked.pl

=head1 SYNOPSIS

  blast_with_ref_cov_forked.pl -- this program will blast all queries in a dir
              against all blastdbs in a dir
    
    for every pair it will also calculate the amnt of non-redundant
    ref sequence covered given a percent identity cutoff.

Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta directory
 --blastdb     Is your blastdb directory
 --config      Is the configuration for your blast
 --outdir      Is your output dir
 --procs       Is the number of forks to run
 --identity    Is the percent identity required for 
                 inclusion in the coverage statistic

ncbi blast executables must be in your path
or blastall must be specified in your config

The config must also specify parameters for blast
except the the query and the reference db

BLASTALL=... (optional)
PARAMETERS=...

=head1 DESCRIPTION

Given a dir of fasta files and a dir of blast dbs, 
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

my ($help, $fasta, $blastdb, $config, $outdir, $procs, $identity);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'b|blastdb=s'     => \$blastdb,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    'i|identity=s'    => \$identity,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $blastdb, $config, $procs, $identity){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/raw`;
`mkdir -p $outdir/parse`;
`mkdir -p $outdir/results`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# get all the queries
opendir (Q, "$fasta");
my @queries = grep (/^.+\..+$/, readdir(Q));
closedir (Q);

# get all the ref genomes
opendir (G, "$blastdb");
my @genomes = grep (/^.+\..+$/, readdir(G));
closedir (G);

my %genomes;
foreach my $file (@genomes){
    my @file = split (/\./, $file);
    pop @file;
    my $store = join (".", @file);
    $genomes{$store} = 1;
}

# note that there is no screen for blast
# against self as this is a needed expt in 
# the next gen sequencing case
my $pm = Parallel::ForkManager->new($procs);
foreach my $query (@queries){
    foreach my $genome (sort keys %genomes){
	$pm->start and next;
	
	# do the forked blast
	warn "Blasting $query $genome\n";
	my $name = blast ($conf, "$fasta/$query", "$blastdb/$genome", "$outdir/raw"); 
	
	# parse the blast result
	warn "Parsing $query $genome\n";
	my ($refdata, $qrydata) = 
	    parse ("$outdir/raw", $name, "$outdir/parse", $identity);

	# tabulate reference covered and queries hit
	warn "Tabulating $query $genome\n";
	my $coverage = genome_coverage ($refdata);
	
	# output
	my $covered ={};
	open (REF, ">$outdir/results/$name.refhits");
	foreach my $chrom (sort keys %$coverage){
	    foreach my $region (sort {$a <=> $b} keys %{$coverage->{$chrom}}){
		my $start = $coverage->{$chrom}->{$region}[0];
		my $end   = $coverage->{$chrom}->{$region}[1];
		my $len   = $end - $start + 1;
		
		$covered->{$chrom} += $len;
		print REF "$chrom\t$region\t$start\t$end\n";
	    }
	}
	close (REF);
	
	open (SUM, ">$outdir/results/$name.refsum");
	foreach my $chrom (sort keys %$covered){
	    print SUM "$chrom\t$covered->{$chrom}\n";
	}
	close (SUM);

	open (QRY, ">$outdir/results/$name.qryhits");
	foreach my $qry (sort keys %$qrydata){
	    foreach my $hit (sort {$a <=> $b} keys %{$qrydata->{$qry}}){
		print QRY "$qry\t$hit\t$qrydata->{$qry}->{$hit}\n";
	    }
	}
	close (QRY);
	$pm->finish;
    }
}
$pm->wait_all_children;

#####SUBS#####

sub genome_coverage {
    my $refdata = shift;
    
    # sort ref data by chrom and by start
    my @refdata = sort {$a->{chrom} cmp $b->{chrom} || 
			    $a->{start} <=> $b->{start}  
		    } @$refdata;
    
    # get the first element
    my $first  = shift @refdata;
    my $chrom1 = $first->{chrom};
    my $start1 = $first->{start};
    my $end1   = $first->{end};

    # do joins on overlapping hits to get the
    # non redundant coverage
    my $coverage = {};
    my $counter  = 0;
    foreach my $point (@refdata){
	
	# get data                                                                                
	my $chrom2 = $point->{chrom};
	my $start2 = $point->{start};
	my $end2   = $point->{end};

	if ($chrom1 eq $chrom2){
	    if (($start1 > $end2) or ($end1 < $start2)) {
		$counter++;

		# they don't overlap at all                                                            
		$coverage->{$chrom1}->{$counter} = [$start1, $end1];
		$chrom1 = $chrom2;
		$start1 = $start2;
		$end1   = $end2;
	    }
	    else {
		my @coor = sort {$a <=> $b} ($start1, $start2, $end1, $end2);
		
		# get the extremes as the new range
		$start1 = shift @coor;
		$end1   = pop @coor;
	    }
	}
	else {
	    $counter = 0; # CHECK WHAT ELSE TO DO HERE
	}

    }
    
    # put in the last values
    $coverage->{$chrom1}->{$counter + 1} = [$start1, $end1];
    
    return ($coverage);

}
    
sub parse {
    my $indir  = shift;
    my $infile = shift;
    my $outdir = shift;
    my $in = Bio::SearchIO->new(-file   => "$indir/$infile",
				-format => 'blast');
    
    open (OUT, ">$outdir/$infile.parse");
    my $qrydata = {};
    my @refdata;
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
#	    last if ($hit_count > 1); # get one hit
	    
	    my $hsp_count = 0;
	    while (my $hsp = $hit->next_hsp ) {
		$hsp_count++;
#		last if ($hsp_count > 1); # get one HSP 
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

		if ($hsp->frac_identical('total') >= $identity){
		    push @refdata, 
		    {
			'chrom' => $hit->name,
			'start' => $hsp->start('hit'),
			'end'   => $hsp->end('hit'),
		    };

		    $qrydata->{$result->query_name}->{$hit_count}++

		}
		
	    }
	    
	}
	
    }
    
    return (\@refdata, $qrydata);
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
