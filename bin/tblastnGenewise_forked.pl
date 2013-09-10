#!/usr/bin/perl -w

=head1 NAME

tblastnGenewise_forked.pl

=head1 SYNOPSIS

  tblastnGenewise_forked.pl -- 
              
Options:

 --help        Show brief help and exit
 --fastadb     Is your query fasta directory (ie, assembly contigs)
    (individual nuc sequences where you want to call genes)
 --aadb        Is your proteins database
    (individual sequences with protein orthologs)
 --config      Is the configuration for your blast and genewise
 --outdir      Is your output dir
 --nblastdb    Is the blastdb of the nuc sequences (ie, assembly contigs)
 --pblastdb    Is the blastdb of the pt sequences (protein ortholog database)
 --procs       Is the number of forks to run

genewise must be in your path.
blast must be in your path.
make sure genewiseconfigdir is defined in your env.

The config must specify parameters for genewise and blast

BLAST=....
GENEWISE=....

=head1 DESCRIPTION

Given a dir of dna seq and a dir of proteins, 
    find genes using tblastn and model with genewise.

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
use Bio::SeqIO;

my ($help, $fastadb, $aadb, $config, $outdir, $procs, $blastdb, $pblastdb, $nblastdb);
GetOptions(
    'h|help'          => \$help,
    'f|fastadb=s'     => \$fastadb,
    'o|outdir=s'      => \$outdir,
    'a|aadb=s'        => \$aadb,
    'c|config=s'      => \$config,
    'p|procs=s'       => \$procs,
    'n|nblastdb=s'    => \$nblastdb,
    'x|pblastdb=s'    => \$pblastdb,
    ) or pod2usage;

pod2usage if $help;

for my $option ($fastadb, $outdir, $aadb, $config, $procs, $nblastdb, $pblastdb){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/tblastn/raw`;
`mkdir -p $outdir/tblastn/parse`;
`mkdir -p $outdir/blastp/raw`;
`mkdir -p $outdir/blastp/parse`;
`mkdir -p $outdir/genewise/raw`;
`mkdir -p $outdir/genewise/proteins`;
`mkdir -p $outdir/genewise/cdnas`;
`mkdir -p $outdir/genewise/gff3`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# get the dna seqs
opendir (F, "$fastadb");
my @fastadb = sort (readdir (F));
shift @fastadb;
shift @fastadb;
closedir (F);

my $fastas = {};
foreach my $fasta (@fastadb){
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fastadb/$fasta");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id    = $sequence_obj->display_id();
	$fastas->{$id} = $fasta;
    }
}

# get the protein sequences
opendir (P, "$aadb");
my @aadb = sort (readdir (P));
shift @aadb;
shift @aadb;
closedir (P);

my $aas = {};
foreach my $aa (@aadb){
    my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$aadb/$aa");
    while (my $sequence_obj = $seqin->next_seq()){
	my $id    = $sequence_obj->display_id();
	$aas->{$id} = $aa;
    }
}


# do the forked gene search using tblastn and genewise
#my $pm = Parallel::ForkManager->new($procs);
foreach my $aa (@aadb){
#    $pm->start and next;
    
    open (PARSE, ">$outdir/$aa.parse");
    
    # run the tblastn and parse
    print STDERR "TBLASTN $aa\n";
    my $nbname = run_blast ($conf, 
			    "TBLASTN", 
			    "$aadb/$aa", 
			    $nblastdb, 
			    "$outdir/tblastn/raw");
    
    my $nbparse = parse_blast ("$outdir/tblastn/raw", 
			       $nbname, 
			       "$outdir/tblastn/parse");
    
    # do genewise searches for every significant blast hit
    print STDERR "GENEWISE $aa\n";
    foreach my $query (sort keys %$nbparse){

	# report and bail of there are no tblastn hits
	if ($nbparse->{$query} eq "NO HITS"){
	    print PARSE "$query\tNo hits\n";
	    next;
	}
	
	foreach my $hit (sort keys %{$nbparse->{$query}}){
	    
	    # run genewise for this hit
	    print STDERR "GW $query $hit\n";
	    my $gname = run_genewise ($conf, 
				      $aadb, 
				      $aas->{$query}, 
				      $fastadb, 
				      $fastas->{$hit}, 
				      "$outdir/genewise/raw");
	    
	    # parse genewise for this hit
	    my $gparse = parse_genewise ("$outdir/genewise/raw", 
					 $gname);
	    
	    # print out gff3 for this hit
	    open (GGFF3, ">$outdir/genewise/gff3/$gname.gff3");
	    my @lines = split (/\n/, $gparse->{'gff3'});
	    foreach my $line (@lines){
		my @line = split (/\t/, $line);
		my $atts = pop (@line);
		my $newatts = "Name=" . $atts . ";" . "Note=" . $gparse->{'data'}[1];
		my $newline = join "\t", @line;
		print GGFF3 "$newline\t$newatts\n";
	    }
	    close (GGFF3);
	    
	    # print out the proteins/cdna and recip blast back to pt reference
	    open (GPT, ">$outdir/genewise/proteins/$gname.fa");
	    print GPT "$gparse->{'protein'}[0]|$gparse->{'data'}[1]\n$gparse->{'protein'}[1]\n";
	    close (GPT);
	    
	    open (GCDNA, ">$outdir/genewise/cdnas/$gname.fa");
	    print GCDNA "$gparse->{'cdna'}[0]|$gparse->{'data'}[1]\n$gparse->{'cdna'}[1]\n";
	    close (GCDNA);
	    
	    print STDERR "BLASTP $query $hit\n";
	    my $pbname = run_blast ($conf, 
				    "BLASTP", 
				    "$outdir/genewise/proteins/$gname.fa", 
				    $pblastdb, 
				    "$outdir/blastp/raw");
	    
	    my $pbparse = parse_blast ("$outdir/blastp/raw", 
				       $pbname, 
				       "$outdir/blastp/parse");

	    # print out the log data
	    # tblastn
	    my $tbnstart = $nbparse->{$query}->{$hit}[0];
	    my $tbnend   = $nbparse->{$query}->{$hit}[1];
	    my $tbnori   = $nbparse->{$query}->{$hit}[2];
	    my $tbneval  = $nbparse->{$query}->{$hit}[3];
	    print PARSE "$query\t$hit\t$tbnstart\t$tbnend\t$tbnori\t$tbneval\t";
	    
	    # genewise
	    my $gwstart  = $gparse->{'data'}[5];
	    my $gwend    = $gparse->{'data'}[6];
	    my $gwbits   = $gparse->{'data'}[0];
	    print PARSE "$gwstart\t$gwend\t$gwbits\t";
	    
	    # recip blastp
	    foreach my $pquery (sort keys %$pbparse){

		# report and bail of there are no tblastn hits                                        
		if ($pbparse->{$pquery} eq "NO HITS"){
		    print PARSE "$pquery\tNo hits\n";
		    next;
		}
		
		my @phits;
		foreach my $phit (sort keys %{$pbparse->{$pquery}}){
		    push (@phits, $phit);
		}
		my $phits = join ";", @phits;
		print PARSE "$pquery\t$phits\n";
		
	    }
	    
	}
    }
    
#    $pm->finish;
}
#$pm->wait_all_children;


#####SUBS#####

sub parse_genewise{
    my $dir = shift;
    my $file = shift;

    local $/ = "//\n";
    
    my $parse = {};
    open (GWO, "$dir/$file");
    while (my $chunk = <GWO>){
	chomp $chunk;
	
	if ($chunk =~m/Bits\s{3}Query/){
	    my ($headers, $data) = split ("\n", $chunk);
	    my @gwdata = split (/\s+/, $data);
	    
	    $parse->{'data'} = [@gwdata];
	}
	elsif ($chunk =~m/\.pep/){
	    my @lines = split (/\n/, $chunk);
	    my $id; 
	    my @seq;
	    foreach my $line (@lines){
		if ($line =~m/^Making/){
		    next;
		}
		elsif ($line =~m/^>/){
		    $id = $line;
		}
		else {
		    push (@seq, $line);
		}
	    }
	    my $seq = join "\n", @seq;
	    $parse->{'protein'} = [$id, $seq];
	}
	elsif ($chunk =~m/\.sp\n/){
	    my @lines =split (/\n/, $chunk);
            my $id = shift (@lines);
            my $seq = join "\n", @lines;
            $parse->{'cdna'} = [$id, $seq];
	}
	elsif ($chunk =~m/genewise-prediction/){
	    $parse->{'gff3'} = $chunk;
	}
	else {
	    next;
	}
    }
    
    return ($parse);
}

sub run_genewise{
    my $conf  = shift;
    my $aadb = shift;
    my $aafile = shift;
    my $fastadb = shift;
    my $fasta = shift;
    my $out   = shift;

    my $cmd = "genewise";
    $cmd .= " $conf->{GENEWISE}";
    $cmd .= " $aadb/$aafile";
    $cmd .= " $fastadb/$fasta";
    $cmd .= " > $out/$aafile.$fasta.out";
    `$cmd`;

    return ("$aafile.$fasta.out");
}

sub parse_blast {
    my $indir  = shift;
    my $infile = shift;
    my $outdir = shift;

    my $in = Bio::SearchIO->new(-file   => "$indir/$infile",
                                -format => 'blast');

    my $parse = {};
    open (F, ">$outdir/$infile.parse");
    while (my $result = $in->next_result){
        my $queryname = $result->query_name;

        # bail if there are no hits                                                                 
        unless ($result->hits){
            print F "$queryname\tNo hits found\n";
            $parse->{$queryname} = 'NO HITS';
            next;
        }
	
	# cycle through the hits
        my @hits = $result->hits;
	foreach my $hit (@hits){
	    my $hitname = $hit->name;

	    # cycle through the hsps and                                                                     
	    # compile and report max and min hsp coordinates
	    my $minstart = 1000000000000000000000000000000;
	    my $maxend   = 0;
	    my $mineval = 10;
	    my $ori;
	    while (my $hsp = $hit->next_hsp){
		my $hspstart = $hsp->start('hit');
		my $hspend   = $hsp->end('hit');
		my $hspori   = $hsp->strand('hit');
		my $hspeval  = $hsp->evalue();
		$hspeval =~s/\,//g;

		($minstart = $hspstart) if ($hspstart < $minstart);
		($maxend   = $hspend) if ($hspend > $maxend);
		($mineval  = $hspeval) if ($hspeval < $mineval);
		$ori = $hspori;
	    }

	    print F "$queryname\t$hitname\t$minstart\t$maxend\t$ori\t$mineval\n";

	    my @data;
	    push (@data, $minstart, $maxend, $ori, $mineval);
	    $parse->{$queryname}->{$hitname} = \@data;
        }
    }
    close (F);
    
    return ($parse);
}

sub run_blast{
    my $conf   = shift;
    my $type   = shift;
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

    my $cmd = "blastall";
    $cmd .= " $conf->{$type}";
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
    close (F);

    return (\%config);
}
