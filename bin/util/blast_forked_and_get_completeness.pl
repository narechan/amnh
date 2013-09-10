#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Bio::SearchIO;

## Stage A: Blast ##

my ($help, $fasta, $blastdb, $config, $outdir, $procs, $delta, $total);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'b|blastdb=s'     => \$blastdb,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
	'd|delta:f'       => \$delta,
	't|total:i'       => \$total
    ) or pod2usage();

pod2usage() if $help;

for my $option ($fasta, $outdir, $blastdb, $config, $procs){
    (warn ("Missing a required option\n") and pod2usage())
	unless ($option);
}


# ...

`mkdir -p $outdir/raw`;
`mkdir -p $outdir/parse`;

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


my @blasted; # for keep on the list of all filename from Blast output


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
    
    # keeps the list of Blasted file. To be checked for his correctness 
    #push @blasted, "$outdir/parse/$fasta.$blastdb.out.parse";
    # Debug
    #print STDERR "Old pushed file was: [$outdir/parse/$name.parse] \n";
    #print STDERR "New pushed file is : [$outdir/parse/$fasta.$blastdb.out.parse] \n";
    # End Debug
    
    
    $pm->finish;
}
$pm->wait_all_children;

# ... 
# read the list of Blasted file that it has found inside the output directory
print STDERR "Read the list of files [$outdir/parse/*.out.parse] \n";
@blasted = grep { -T } glob "$outdir/parse/*.out.parse";

# loop on all files, made by Blast
for my $item (@blasted) {
	print STDERR "Started the working on the file $item\n";
	my $records_ref = &stage_B($item);
	
	my $delta_default = 0.8;
	my $total_default = scalar @{$records_ref};

	$delta = $delta || $delta_default;
	$total = $total || $total_default;

	&stage_C($records_ref, $delta, $total);
}


##########################################
#               S u b s                  #
##########################################

sub stage_B() {
	my $in1 = shift;
	
	open IN1, $in1 or die("$!\n");
	my @rows = grep { !/^\w+\s+\d+\s+No\s+/ } <IN1> ;
	close IN1;

	my @rr =  map { join "\t", (split)[0,2,1,3,7] } @rows;

	@rr = sort { 
				(split /\t/, $a)[4] <=>
				(split /\t/, $b)[4] 
			} @rr;

	my %rr = map { (split /\t/, $_)[0] => $_ } @rr;

	#print map { $rr{$_} ,"\n" } sort keys %rr;

	@rr = map { $rr{$_} } sort keys %rr; 

	# debug
	print map { $_, "\n"} @rr;
	
	return \@rr;
}

sub stage_C() {
	my @rr = @{shift()};
	my $delta = shift;
	my $total = shift;

	die "Error: Total $total should be > 0" if defined $total && $total <= 0;
							  
	my $line;
	my %hits;
	my $n = 0;

	for $line (@rr) {
		my @columns = split(/\t/, $line);
		my $num_columns = scalar(@columns);
	   
		my ( $query, $subject, $query_len, $subject_len, $bitscore ) = split(/\t/, $line);

		# assumes subject lengths are valid integers
		# zero subject lengths are invalid, i.e., are not hits
		next unless $subject_len > 0;

		$n++;

		# check each input line for 5 consistent elements
		die "Error: expected 5 elements in line $n, but found $num_columns\n" if $num_columns != 5;

		my $coverage = $query_len / $subject_len;

		if ( ! defined $hits{ $subject }{'bitscore'} || $bitscore > $hits{ $subject }{'bitscore'} )
		{
			$hits{ $subject }{'bitscore'} = $bitscore;
			$hits{ $subject }{'coverage'} = $coverage;
		}

		
	}
	
	my $count = 0;
	for my $subject (keys %hits) { $count++ if $hits{ $subject }{'coverage'} >= $delta };

	my $running = scalar(keys %hits);
	my $N = defined $total ? $total : $running;

	# make sure the total $N is not less than the running total $running
	die "Error: did not expect running total $running to be greater than total $total\n" if ($running > $N);

	my $local = 100 * ( $count / $running  );
	my $completeness = 100 * ( $count / $N );


	# Using a Format  so the output will be more readable
	$n = sprintf("%d", $n);
	$N = sprintf("%d",$N);
	$running = sprintf("%d",$running);
	$count = sprintf("%d",$count);
	$local = sprintf("%.1f",$local);
	$completeness = sprintf("%.1f",$completeness);


format STDOUT =
+==========================================================================================================+
| Entries        Total          Running Total   Number > @<<<<<    Local Completeness   Total Completeness |
$delta
+----------------------------------------------------------------------------------------------------------+
| @<<<<<<<<<<<<  @<<<<<<<<<<<<  @<<<<<<<<<<<<   @<<<<<<<<<<<<      @<<<<<<<<<<<<        @<<<<<<<<<<<<      |
$n, $N, $running, $count, $local, $completeness
+==========================================================================================================+

.

	write;

}


##########################################
#  Subs from blast_forked.pl unchanged   #
##########################################

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




__END__

=head1 NAME

blast_forked_and_get_completeness.pl

=head1 SYNOPSIS

  blast_forked_and_get_completeness.pl -- 
              

Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta directory or fasta file (script can handle both). Required
 --blastdb     Is your blastdb. Required
 --config      Is the configuration for your blast. Required
 --outdir      Is your output dir. Required
 --procs       Is the number of forks to run. Required
 --delta       Is the coverage percent cutoff (0.0-1.0) (default 0.8). Optional
 --total       Is the total number of reference transcripts (>1) (default the number of rows). Optional
 
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

Modified by Leo Manfredi, manfredi@cpan.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

