#!/usr/bin/perl -w

=head1 NAME

nucmer_forked.pl

=head1 SYNOPSIS

  nucmer_forked.pl -- this program will globally align (using mummer)
              a set of query genomes against a set of ref genomes

Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta directory
 --ref         Is your ref genome directory
 --config      Is the configuration for your blast
 --outdir      Is your output dir
 --procs       Is the number of forks to run

mummer executables must be in your path

The config must also specify parameters for mummer executables

NUCMER=...
SHOW-COORDS=...

=head1 DESCRIPTION

Given a dir of query genomes and a dir of ref genomes
    do and parse the mummer alignments

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

my ($help, $fasta, $ref, $config, $outdir, $procs);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'r|ref=s'         => \$ref,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $ref, $config, $procs){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/raw`;
`mkdir -p $outdir/results`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# get all the queries
opendir (Q, "$fasta");
my @queries = sort (readdir(Q));
shift @queries;
shift @queries;
closedir (Q);

# get all the ref genomes
opendir (G, "$ref");
my @genomes = sort (readdir(G));
shift @genomes;
shift @genomes;
closedir (G);

# do the forked mummer algns
my $pm = Parallel::ForkManager->new($procs);
foreach my $query (@queries){
    foreach my $genome (@genomes){
	$pm->start and next;
	
	warn "Aligning $query $genome\n";
	mummer ($conf, "$fasta/$query", "$ref/$genome", "$outdir/raw/$query.$genome"); 
	
	warn "Parsing $query $genome\n";
	my $refdata = parse ("$outdir/raw/$query.$genome");
	my $coverage = genome_coverage ($refdata);
	
	# output                                                                                   
	my $covered ={};
        open (REF, ">$outdir/results/$query.$genome.refhits");
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

        open (SUM, ">$outdir/results/$query.$genome.refsum");
        foreach my $chrom (sort keys %$covered){
            print SUM "$chrom\t$covered->{$chrom}\n";
	    print STDERR "RefCov\t$query\t$genome\t$covered->{$chrom}\n";
        }
        close (SUM);

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


sub parse{
    my $file = shift;
   
    my @refdata;
    open (F, "$file.showcoords");
    while (my $line = <F>){
	chomp $line;
	my @data = split (/\t/, $line);
	
	push @refdata,
	{
#	    'chrom' => $data[11],
	    'chrom' => $data[7],
	    'start' => $data[0],
	    'end'   => $data[1],
	};
    }
    close (F);
    
    return (\@refdata);
}
	

sub mummer{
    my $conf   = shift;
    my $fasta  = shift;
    my $ref    = shift;
    my $out    = shift;

    my $fasta_name;
    if ($fasta =~/\//g){
	$fasta =~m/.*\/(.*)$/;
	$fasta_name = $1;
    }

    else {
	$fasta_name = $fasta;
    }

    my $ref_name;
    if ($ref =~/\//g){
	$ref =~m/.*\/(.*)$/;
	$ref_name = $1;
    }

    else {
	$ref_name = $ref;
    }

    # nucmer
    my $nucmer = "nucmer";
    ($nucmer .= " $conf->{'NUCMER'}") if ($conf->{'NUCMER'});
    $nucmer .= " -p $out";
    $nucmer .= " $ref";
    $nucmer .= " $fasta";
    `$nucmer`;
    
    # show-coords
    my $showcoords = "show-coords";
    ($showcoords .= " $conf->{'SHOW-COORDS'}") if ($conf->{'SHOW-COORDS'});
    $showcoords .= " $out.delta";
    $showcoords .= " > $out.showcoords";
    `$showcoords`;

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
