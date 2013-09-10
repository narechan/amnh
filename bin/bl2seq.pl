#!/usr/bin/perl -w

=head1 NAME

bl2seq.pl

=head1 SYNOPSIS

  bl2seq -- 
              
Options:

 --help        Show brief help and exit
 --fasta       Is your query fasta file
 --config      Is the configuration for your blast
 --outdir      Is your output dir

ncbi blast executables must be in your path
or bl2seq must be specified in your config

The config must also specify parameters for blast
except the the query and the reference seqs

BLASTALL=... (optional)
PARAMETERS=...

=head1 DESCRIPTION

Given a fasta file blast all members in pairwise fashion
and parse.

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
use Bio::SeqIO;
use Bio::AlignIO;

my ($help, $fasta, $config, $outdir);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fasta, $outdir, $config){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/raw`;
`mkdir -p $outdir/parse`;
`mkdir -p $outdir/alns`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# parse the fasta file
my $sequence = {};
my $sin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fasta");
while (my $sequence_obj = $sin->next_seq()){
    my $id  = $sequence_obj->display_id();
    my $seq = $sequence_obj->seq();
    $sequence->{$id} = $seq;
}

# do the blast
foreach my $id1 (sort (keys %$sequence)){
    foreach my $id2 (sort (keys %$sequence)){
	next if ($id1 eq $id2);

	# parse deflines and only blast 
	# those seqs in same trinity comp
	# comment out if want to do all pw
	my @id1 = split (/\_/, $id1, 2);
	my @id2 = split (/\_/, $id2, 2);
	my $one = shift @id1;
	my $two = shift @id2;
	next unless ($one eq $two);
	
	print STDERR "Working on $id1 $id2\n";
	
	# generate tmp fasta files
	open (ONE, ">/tmp/$id1.fa");
	open (TWO, ">/tmp/$id2.fa");
	print ONE ">$id1\n$sequence->{$id1}\n";
	print TWO ">$id2\n$sequence->{$id2}\n";
	
	# bl2seq
	my $filename   = blast ($conf, "/tmp/$id1.fa", "/tmp/$id2.fa", "$outdir/raw");
	my $hspstokeep = parse ("$outdir/raw", $filename, "$outdir/parse", $id1, $id2);

	# get the blast hsp alignments
	my $counter = 0;
	foreach my $hsp (@$hspstokeep){
	    $counter++;
	    
	    my $aln = $hsp->get_aln;
	    my $alnIO = Bio::AlignIO->new(-format =>"fasta", -file => ">$outdir/alns/$filename.$counter.fa");
	    $alnIO->write_aln($aln);
	}

	# rm the tmp files
	`rm /tmp/$id1.fa`;
	`rm /tmp/$id2.fa`;
    }
    
    # don't need both directions so pare
    # off the top
    delete $sequence->{$id1};
}
    
#####SUBS#####

sub parse {
    my $indir  = shift;
    my $infile = shift;
    my $outdir = shift;
    my $id1    = shift;
    my $id2    = shift;
    my $in = Bio::SearchIO->new(-file   => "$indir/$infile",
				-format => 'blast');

    my @hspstokeep;
    open (OUT, ">$outdir/$infile.parse");
    while( my $result = $in->next_result ) {
        my $hit_count = 0;  
        unless ($result->hits) {
            print OUT join(
			   "\t",
			   $result->query_name,        #1
			   $result->query_length,      #2
			   $id2,                       #2.5 (label with hit name even if not sig)
			   'No significant hits found',            #3
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
		
		push (@hspstokeep, $hsp);
		
		print OUT join(
			       "\t",
			       $result->query_name,           #1
			       $result->query_length,         #2
			       $hit->name,                    #2.5 (repeat to keep columns)
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
    return (\@hspstokeep);
}

sub blast{
    my $conf   = shift;
    my $onefile  = shift;
    my $twofile  = shift;
    my $out    = shift;

    my $onefile_name;
    if ($onefile =~/\//g){
	$onefile =~m/.*\/(.*)$/;
	$onefile_name = $1;
    }

    else {
	$onefile_name = $onefile;
    }

    my $twofile_name;
    if ($twofile =~/\//g){
        $twofile =~m/.*\/(.*)$/;
        $twofile_name = $1;
    }

    else {
        $twofile_name = $twofile;
    }

    my $cmd = $conf->{BL2SEQ};
    $cmd .= " $conf->{PARAMETERS}";
    $cmd .= " -i $onefile";
    $cmd .= " -j $twofile";
    $cmd .= " -o $out/$onefile_name.$twofile_name.out";
    `$cmd`;
    
    return ("$onefile_name.$twofile_name.out");
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
    if (! exists($config{"BL2SEQ"})){
        my $blastall = `which bl2seq`;
        chomp $blastall;

        if (defined ($blastall)){
            $config{"BL2SEQ"} = $blastall;
        }
        else {
            die "no BL2SEQ specified!\n";
	}
    }

    return (\%config);
    close (F);
}
