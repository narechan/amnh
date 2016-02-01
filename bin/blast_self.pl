#!/usr/bin/perl -w

=head1 NAME

blast_self.pl

=head1 SYNOPSIS

  blast_self.pl -- 
              

Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta file
 --config      Is the configuration for your blast
 --identity    Is the percent identity above which regions are screened   
 --outdir      Is your output dir
 --procs       Is the number of forks to run

ncbi blast executables must be in your path
or blastall must be specified in your config

The config must also specify parameters for blast
except the the query and the reference db

PARAMETERS=...

=head1 DESCRIPTION

Given a fasta file, create a blast database of that file
    and blast against self to generate masked regions

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

my ($help, $fasta, $config, $outdir, $procs, $ident);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    'i|identity=s'    => \$ident,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $config, $procs, $ident){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/raw`;
`mkdir -p $outdir/parse`;
`mkdir -p $outdir/substrings`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# build the blast database
`formatdb -i $fasta -o T -p F -V`;

# slice up the sequence to create the queries
my $txt;
open (F, "$fasta");
while (my $line = <F>){
    chomp $line;
    next if ($line =~m/^>/);
    $txt .= $line;
}
my @substrings;
my $start = 0;
while (length $txt > 0){
    $start++;
    my $end = $start + 100 - 1;

    open (S, ">$outdir/substrings/$start-$end.fa");
    my $substring = substr ($txt, 0, 100, '');
    print S ">$start-$end\n$substring\n";
    close (S);

    push (@substrings, "$start-$end.fa");
    $start = $end;
}

# do the forked blast
my $pm = Parallel::ForkManager->new($procs);
foreach my $query (@substrings){
    $pm->start and next;
    
    my ($range, $crap) = split (/\./, $query);
    print STDERR "Blasting $query\n";
    my $name = blast ($conf, "$outdir/substrings/$query", "$fasta", "$outdir/raw"); 
    
    $pm->finish;
}        
$pm->wait_all_children;

# do the parse
opendir (D, "$outdir/raw");
my @raws = sort (readdir (D));
shift @raws;
shift @raws;
closedir (D);

my $regions = {};
foreach my $name (@raws){
    print STDERR "Parsing $name\n";
    my @name = split (/\./, $name);
    my ($qstart, $qend) = split (/\-/, $name[0]);
    
    my $in = Bio::SearchIO->new(-file   => "$outdir/raw/$name",
				-format => 'blast');
    my $result = $in->next_result;
    my $hit_count = 0;
    while( my $hit = $result->next_hit ) {
	$hit_count++;
	
	my $hsp_count = 0;
	while (my $hsp = $hit->next_hsp ) {
	    $hsp_count++;
#	    last if ($hsp_count > 1); # get one HSP
	    
	    my $frac   = $hsp->frac_identical;
	    my $hstart = $hsp->start('hit');
	    my $hend   = $hsp->end('hit');
	    
	    if (($hstart >= $qstart) and ($hend <= $qend)){
		print STDERR "$hit_count\t$hsp_count\tself\n";
		next;
	    }
	    if ($frac >= $ident){
		print STDERR "$hit_count\t$hsp_count\tstored\n";
		$regions->{$name[0]} = 1;
	    }
	}
    }
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
#	    last if ($hit_count > 1); # get one hit
	    
	    my $hsp_count = 0;
	    while (my $hsp = $hit->next_hsp ) {
		$hsp_count++;
#		last if ($hsp_count > 1); # get one HSP 

		# length requirement if desired
#		my $hsplen = $hsp->hsp_length;
#		my $qrylen = $result->query_length;
#		my $lenreq = $hsplen / $qrylen;
#		if (($lenreq < 0.95) or ($lenreq > 1.05)){
#		    print OUT join(
#				   "\t",
#				   $result->query_name,        #1                                        
#				   $result->query_length,      #2                                            
#				   'Length insufficient',            #3                             
#				   '0',                        #4                                            
#				   '1',                        #5                                              
#				   '1',                        #6                                             
#				   '1000',                     #7                                              
#				   '0',                        #8                                           
#				   '0',                        #9                                            
#				   '0',                        #10                                           
#				   '0',                        #11                                           
#				   '0',                        #12                                        
#				   '0',                        #13                                           
#				   '0',                        #14                                           
#				   '0',                        #15                                           
#				   '0',                        #16                                             
#				   '0',                        #17                                         
#				   '0',                        #18                                         
#				   '0',                        #19                                          
#				   $result->query_description, #20                                           
#				   'NULL',                     #21                                            
#				   ), "\n";
#		    next;
#		}
		    
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
