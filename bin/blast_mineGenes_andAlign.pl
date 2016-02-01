#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SearchIO;
use Bio::Index::Fasta;

my ($help, $fasta, $qrydir, $config, $gff3dir, $genomefadir, $outdir);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta, #this is the pt you are blasting against the db
    'o|outdir=s'      => \$outdir,
    'q|qrydir=s'       => \$qrydir, #this is a directory of proteomes to query against
    'c|config=s'      => \$config,
    'g|gffdir=s'         => \$gff3dir, #this is a directory of gffs corresponding to those proteomes
    'x|genomefadir=s'    => \$genomefadir, #this is a directory of genomes corresponding to those proteomes
    ) or pod2usage;

# note that all file names in the three dirs above need to be the same (even suffixes), and if there are uneven
# numbers of files, the genomefadir should contain the most complete set (usually the case anyway)

pod2usage if $help;

for my $option ($fasta, $outdir, $qrydir, $config, $gff3dir, $genomefadir){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/blastp`;
`mkdir -p $outdir/tblastn`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

opendir (D, "$genomefadir");
my @qryfas = sort (readdir (D));
shift @qryfas;
shift @qryfas;
closedir (D);

foreach my $qryfa (@qryfas){
    `cp $qrydir/$qryfa $outdir/$qryfa.proteome`;
    `cp $genomefadir/$qryfa $outdir/$qryfa.genome`;
    
    # index the pt qry fasta
    my $index = Bio::Index::Fasta->new(-filename => $qryfa . ".idx", -write_flag => 1);
    $index->make_index("$outdir/$qryfa.proteome");

    # index the genome fasta
    my $gindex = Bio::Index::Fasta->new(-filename => $qryfa . ".idx", -write_flag => 1);
    $gindex->make_index("$outdir/$qryfa.genome");

    # format the pt qry fasta and nuc query fa and move stuff to local
    `formatdb -i $outdir/$qryfa.proteome -o T -p T -V`;
    `formatdb -i $outdir/$qryfa.genome -o T -p F -V`;

    # parse the gff file and key on pt names with genomic locations
    my $gff3coords = {};
    open (G, "$gff3dir/$qryfa");
    my $counter = 0;
    while (my $line = <G>){
	chomp $line;
	$counter++;
	
	next if ($line =~m/^\#/);
	
	my ($chrom,
	    $source,
	    $type,
	    $start,
	    $end,
	    $score,
	    $strand,
	    $phase,
	    $attr) = split (/\t/, $line);
	
	my $tstrand;
	if ($strand eq "+"){
	    $tstrand = 1;
	}
	else {
	    $tstrand = -1;
	}
	
	
	# bail if the type is not cds                                                                          
	next unless ($type eq "CDS");
	
	my @attrs = split (/\;/, $attr);
	foreach my $att (@attrs){
	    my ($key, $value) = split (/\=/, $att);
	    if ($key eq "Name"){
		my $coordstring = $chrom . ":" . $tstrand . ":" . $start . ":" . $end;
		$gff3coords->{$value} = $coordstring;
	    }
	    else {
		next;
	    }
	}
    }

    # do the pt blast to the query proteome
    my $name = blast ($conf, "$fasta", "$outdir/blastp/$qryfa", "$outdir", "BLASTP"); 

    # parse out the best hit and output the sequence
    my $in = Bio::SearchIO->new(-file   => "$outdir/$name",
				-format => 'blast');
    while( my $result = $in->next_result ){
	if ($result->hits){
	    while(my $hit = $result->next_hit){

		# get hit attributes
		my $hname = $hit->name;
		my $hsig  = $hit->significance;
		
		# get pt sequence
		my $sequence = $index->fetch($hname);
		my $id       = $sequence->id();
		my $seq      = $sequence->seq();
		
		open (P, ">$outdir/blastp.pt.seq");
		print P ">$qryfa-$hname\n$seq\n";
		close (P);
	    
		# get nuc sequence
		my ($contig, $strand, $start, $end) = split (/\:/, $gff3coords->{$id});
		my $location = Bio::Location::Simple->new(-start  => $start,
							  -end    => $end,
							  -strand => $strand);
		my $gsequence = $gindex->fetch($contig);
		my $gseq      = $gsequence->subseq($location);
	    
		open (G, ">$outdir/blastp.nuc.seq");
		print G ">$qryfa-$hname\n$gseq\n";
		close (G);
	    
		open (S, ">$outdir/blastp.annot.seq");
		print S "$qryfa\t$hname\t$hsig\t";
		while (my $hsp = $hit->next_hsp){
		    my $hsplen = $hsp->hsp_length;
		    my $numgaps = $hsp->gaps;
		    my $numident = $hsp->num_identical;
		    my $nucdiv = ($hsplen - $numident) / ($hsplen - $numgaps);
		    print S "$nucdiv\n";
		    last;
		}
		close (S);
		last;
	    }
	}
	else{
	    print STDERR "$qryfa\tNo hits\n";
	}
    }

    # do the the tblastn directly to the query contigs
    my $tname = blast ($conf, "$fasta", "$outdir/tblastn/$qryfa", "$outdir", "TBLASTN");

    # parse out the best hit and output the sequence                                                      
    my $tin = Bio::SearchIO->new(-file   => "$outdir/$tname",
				 -format => 'blast');
    while( my $result = $tin->next_result ){
	if ($result->hits){
	    while(my $hit = $result->next_hit){
		
		# get hit attributes                                                                         
		my $hname = $hit->name;
		my $hsig  = $hit->significance;
		
		# get nuc sequence                                                                            
		open (T, ">$outdir/tblastn.annot.seq");
		print T "$qryfa\t$hname\t$hsig\t";
		while (my $hsp = $hit->next_hsp){
		    my $hspst = $hsp->start('hit');
		    my $hspend = $hsp->end('hit');
		    my $hspstrand = $hsp->strand('hit');
		    my $location = Bio::Location::Simple->new(-start  => $hspst,
							      -end    => $hspend,
							      -strand => $hspstrand);
		    my $tsequence = $gindex->fetch($hname);
		    my $tseq      = $tsequence->subseq($location);
		    open (TG, ">$outdir/tblastn.nuc.seq");
		    print TG ">$qryfa-$hname\n$tseq\n";
		    close (TG);
		    
		    my $hsplen = $hsp->hsp_length;
		    my $numgaps = $hsp->gaps;
		    my $numident = $hsp->num_identical;
		    my $nucdiv = ($hsplen - $numident) / ($hsplen - $numgaps);
		    print T "$nucdiv\n";
		    last;
		}
		close (T);
		last;
	    }
	}
	else{
	    print STDERR "$qryfa\tNo hits\n";
	}
    }
    
#`rm -rf $outdir`;
}	    
#####SUBS#####

sub blast{
    my $conf   = shift;
    my $fasta  = shift;
    my $db     = shift;
    my $out    = shift;
    my $prog   = shift;

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
    $cmd .= " $conf->{$prog}";
    $cmd .= " -i $fasta";
    $cmd .= " -d $db";
    $cmd .= " -o $out/$fasta_name.$db_name.$prog.out";
    `$cmd`;
    
    return ("$fasta_name.$db_name.$prog.out");
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
