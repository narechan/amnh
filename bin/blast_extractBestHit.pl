#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SearchIO;
use Bio::Index::Fasta;

#my ($help, $fasta, $qryfa, $config, $gff3, $genomefa, $outdir);
my ($help, $fasta, $config, $genomefa, $outdir);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta, #this is the pt you are blasting against the db
    'o|outdir=s'      => \$outdir,
#    'q|qryfa=s'       => \$qryfa, #this is the pt db
    'c|config=s'      => \$config,
#    'g|gff=s'         => \$gff3, #this is the genome annotation
    'x|genomefa=s'    => \$genomefa, #this is the genome
    ) or pod2usage;

pod2usage if $help;

#for my $option ($fasta, $outdir, $qryfa, $config, $gff3, $genomefa){
for my $option ($fasta, $outdir, $config, $genomefa){ 
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# index the pt qry fasta
#my $index = Bio::Index::Fasta->new(-filename => $qryfa . ".idx", -write_flag => 1);
#$index->make_index($qryfa);

# index the genome fasta
`cp $genomefa $outdir`;
my $genomefa_name;
if ($genomefa =~/\//g){
    $genomefa =~m/.*\/(.*)$/;
    $genomefa_name = $1;
}
else {
    $genomefa_name = $genomefa;
}
my $gindex = Bio::Index::Fasta->new(-filename => "$outdir/$genomefa_name" . ".idx", -write_flag => 1);
$gindex->make_index("$outdir/$genomefa_name");

# format the pt qry fasta and nuc query fa and move stuff to local
#`formatdb -i $qryfa -o T -p T -V`;
`formatdb -i $outdir/$genomefa_name -o T -p F -V`;
#`mv $qryfa.* $outdir`;
#`mv $genomefa.* $outdir`;
=head
# parse the gff file and key on pt names with genomic locations
my $gff3coords = {};
open (G, "$gff3");
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
=cut
=head
# get names
my $qryfa_name;
if ($qryfa =~/\//g){
    $qryfa =~m/.*\/(.*)$/;
    $qryfa_name = $1;
}
else {
    $qryfa_name = $qryfa;
}
=cut
=head
# do the pt blast to the query proteome
my $name = blast ($conf, "$fasta", "$outdir/$qryfa_name", "$outdir", "BLASTP"); 

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
	    print P ">$qryfa_name-$hname\n$seq\n";
	    close (P);
	    
	    # get nuc sequence
	    my ($contig, $strand, $start, $end) = split (/\:/, $gff3coords->{$id});
	    my $location = Bio::Location::Simple->new(-start  => $start,
						      -end    => $end,
						      -strand => $strand);
	    my $gsequence = $gindex->fetch($contig);
	    my $gseq      = $gsequence->subseq($location);
	    
	    open (G, ">$outdir/blastp.nuc.seq");
	    print G ">$qryfa_name-$hname\n$gseq\n";
	    close (G);
	    
	    open (S, ">$outdir/blastp.annot.seq");
	    print S "$qryfa_name\t$hname\t$hsig\t";
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
	print STDERR "$qryfa_name\tNo hits\n";
    }
}
=cut
# do the the tblastn directly to the query contigs
#my $tname = blast ($conf, "$fasta", "$outdir/$genomefa_name", "$outdir", "TBLASTN");
#my $tname = blast ($conf, "$fasta", "$outdir/$genomefa_name", "$outdir", "TBLASTX");
my $tname = blast ($conf, "$fasta", "$outdir/$genomefa_name", "$outdir", "BLASTN"); 

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
#	    open (T, ">$outdir/tblastn.annot.seq");
#	    open (T, ">$outdir/tblastx.annot.seq"); 
	    open (T, ">$outdir/blastn.annot.seq");
            print T "$genomefa_name\t$hname\t$hsig\t";
            while (my $hsp = $hit->next_hsp){
		my $hspst = $hsp->start('hit');
		my $hspend = $hsp->end('hit');
		my $hspstrand = $hsp->strand('hit');
		my $location = Bio::Location::Simple->new(-start  => $hspst,
							  -end    => $hspend,
							  -strand => $hspstrand);
		my $tsequence = $gindex->fetch($hname);
		my $tseq      = $tsequence->subseq($location);
#		open (TG, ">$outdir/tblastx.nuc.seq");
		open (TG, ">$outdir/blastn.nuc.seq");
		print TG ">$genomefa_name-$hname\n$tseq\n";
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
        print STDERR "$genomefa_name\tNo hits\n";
    }
}

#`rm -rf $outdir`;
	    
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
