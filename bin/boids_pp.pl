#!/usr/bin/perl
=head1 NAME                                                                       
                                                                                  
boids_pp.pl
                                                                                  
=head1 SYNOPSIS                                                                   
                                                                                  
  boids_pp.pl --                                                           
    Post processes a set of bootstrapped boids simulations to
    produce robust gene clusters.
                                                                                  
Options:                                                                          
                                                                                  
 --help        Show brief help and exit                                           
 --matrix      Is your nexus file of sequences (boids are partitions)
 --infile      Is the set of boids coordinates you want to analyze
 --xrange      Is the xrange coordinates to cut out; comma separated, no spaces
 --yrange      Is the yrange coordinates to cut out; comma separated, no spaces   
 --hmms        Is a directory of hmms for each of the partitions in the matrix
    (optional; if not provided an hmm library will be created
 --refgb       Is your reference genbank file
 --reffa       Is your reference fasta file
 --outdir      Is your output dir                                                 
 --conf        Is your configuration for embedded programs (hmmbuild,..)
 --score       Is the hmm bitscore cutoff for assignment

=head1 DESCRIPTION                                                                
                                                                                  
Parse flocks a derive clusters
                                                                                  
=head1 AUTHOR                                                                     
                                                                                  
Apurva Narechania                                                                 
anarechania *a|t* amnh.org                                                        
                                                                                  
=head1 COPYRIGHT                                                                  
                                                                                  
Copyright (c) 2012 American Museum of Natural History                             
                                                                                  
This library is free software;                                                    
you can redistribute it and/or modify                                             
it under the same terms as Perl itself.                                           
                                                                                  
=cut                                                                              

# ----------------------------------------------------                            

#####                                                                             
#####                                                                             

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Index::Fasta;
use Ild;

my ($help, $matrix, $outdir, $infile, $xrange, $yrange, $hmmdir, $refgb, $reffa, $configfile, $score);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrix,
    'o|outdir=s'      => \$outdir,
    'i|infile=s'      => \$infile,
    'x|xrange=s'      => \$xrange,
    'y|yrange=s'      => \$yrange,
    'h|hmms=s'        => \$hmmdir,
    'g|refgb=s'       => \$refgb,
    'f|reffa=s'       => \$reffa,	   
    'c|conf=s'        => \$configfile,
    's|score=s'       => \$score,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($matrix, $outdir, $infile, $xrange, $yrange, $refgb, $reffa, $configfile, $score){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir/hmms`;
`mkdir -p $outdir/hmmsearch`;
`mkdir -p $outdir/circos`;

####MAIN####

# parse the config file
my $config = parse_config ($configfile);

# parse the nexus file and store the data
print STDERR "Reading matrix\n";
my $ildobj = Ild->new;
$ildobj->load_aln ($matrix);
my $charsets = $ildobj->get_charsets;
my ($partitions, $lengths) = $ildobj->store_alignment_seq($matrix, $charsets);

# read in hmms or create the hmm library from the nexus file
print STDERR "Processing HMMs\n";
my @hmms;
if ($hmmdir){
    opendir (H, "$hmmdir");
    @hmms = sort (readdir (H));
    shift @hmms;
    shift @hmms;
    closedir (H);
}
else {
    $hmmdir = "$outdir/hmms";
    foreach my $part (sort keys %$partitions){
	print STDERR "Building $part HMM\n";
	open (F, ">$outdir/$part.fa");
	foreach my $taxon (sort keys %{$partitions->{$part}}){
	    my $seq = $partitions->{$part}->{$taxon};
	    $seq =~s/\?/\-/g;
	    print F ">$taxon\n$seq\n";
	}
	close (F);
	
	my $hmmbuild = "hmmbuild";
	$hmmbuild   .= " $config->{'HMMBUILD'}";
	$hmmbuild   .= " $outdir/hmms/$part.hmm $outdir/$part.fa";
	`$hmmbuild`;
	push (@hmms, "$part.hmm");
	`rm $outdir/$part.fa`;
    }
}

print STDERR "Processing the reference\n";    
# index the fasta file                                                              
my $index = Bio::Index::Fasta->new(-filename => $reffa . ".idx", -write_flag => 1);
$index->make_index($reffa);

# parse the genbank file
my $ref = {};
my $seqin = Bio::SeqIO->new(-file=>$refgb, -format=>'Genbank');
my $seqobj = $seqin->next_seq;
my $chrom  = $seqobj->display_id;

# cycle through the features                                                 
open (P, ">$outdir/pt.fa");
#open (N, ">$outdir/nuc.fa");
foreach my $feat ($seqobj->get_SeqFeatures){
    my $type = $feat->primary_tag;
    
    # configure the circos karyotype
    if ($type eq "source"){
	open (K, ">$outdir/circos/karyotype.data");
	my $s = $feat->start;
	my $e = $feat->end;
	print K "chr - $chrom $chrom $s $e gpos25\n";
	close (K);
    }
    
    next unless ($feat->primary_tag =~m/CDS|tRNA|rRNA|transcript/);

    # parse the data
    my $name;
    if ($feat->has_tag("gene")){
	my @names = $feat->get_tag_values("gene");
	$name = join "_", @names;
    }
    elsif ($feat->has_tag("locus_tag")){
	my @names = $feat->get_tag_values("locus_tag");
	$name = join "_", @names;
    }
    else {
	print STDERR "No name!\n";
	die;
    }
    
    my $desc;
    if ($feat->has_tag("product")){
	my @names = $feat->get_tag_values("product");
	$desc = join "_", @names;
    }
    else {
	$desc = "NODESC";
    }
    
    my $start = $feat->start;
    my $end   = $feat->end;
    my $strand = $feat->strand;
    
    # harvest sequence
    my $location = Bio::Location::Simple->new(-start  => $start,
					      -end    => $end,
					      -strand => $strand);
    my $sequence = $index->fetch($chrom);
    my $subseq   = $sequence->subseq($location);

    my $subseq_obj = Bio::Seq->new(-seq => $subseq,
				   -alphabet => 'dna');
    my $ptsubseq_obj = $subseq_obj->translate(-frame=>0);
    my $pt = $ptsubseq_obj->seq;
    print P ">$name\n$pt\n";
#    print N ">$name\n$subseq\n";

    # store data
    $ref->{$name}->{'start'} = $start;
    $ref->{$name}->{'end'} = $end;
    $ref->{$name}->{'strand'} = $strand;
    $ref->{$name}->{'dna'} = $subseq;
    $ref->{$name}->{'protein'} = $pt;
    $ref->{$name}->{'desc'} = $desc;
}

# parse the boids coordinates file and find relavent groups
print STDERR "Finding query groups\n";
my $boidskept = {};
my ($x1, $x2) = split (/\,/, $xrange);
my ($y1, $y2) = split (/\,/, $yrange);
open (B, "$infile");
while (my $line = <B>){
    chomp $line;
    my ($gene, $x, $y) = split (/\t/, $line);
    if (($x >= $x1) and ($x <= $x2) and ($y >= $y1) and ($y <= $y2)){
	$boidskept->{$gene} = 1;
#	print "$gene\n";
    }
    else {
	next;
    }
}
close (B);

# score the groups harvested against the reference using hmmsearch
#print STDERR "Scoring group HMMs against refereence proteins\n";
print STDERR "Scoring group HMMs against refereence genes\n";
open (O, ">$outdir/correspondences.out");
open (ANNOT, ">$outdir/circos/features.annot");
foreach my $boid (sort keys %$boidskept){
    print STDERR "Scoring $boid\n";
    my $hmmsearch = "hmmsearch";
    $hmmsearch   .= " $config->{'HMMSEARCH'}";
    $hmmsearch   .= " -o $outdir/hmmsearch/$boid.hmm.out --tblout $outdir/hmmsearch/$boid.hmm.tblout $hmmdir/$boid.hmm $outdir/pt.fa";
#    $hmmsearch   .= " -o $outdir/hmmsearch/$boid.hmm.out --tblout $outdir/hmmsearch/$boid.hmm.tblout $hmmdir/$boid.hmm $outdir/nuc.fa";
    `$hmmsearch`;
    
    my $refseq;
    my $refscore;
    open (A, "$outdir/hmmsearch/$boid.hmm.tblout");
    while (my $line = <A>){
	chomp $line;
	next if ($line =~m/\#/);
	my @line = split (/\s+/, $line);
	
	if ($line[5] > $score){
	    $refseq = $line[0];
	    $refscore = $line[5];
	}
	else {
	    $refseq = "NONE_SIG";
	    $refscore ="NA";
	}
	last; #just want the top hit
    }
    close (A);
    
    if ($refseq){
	if ($refseq eq "NONE_SIG"){
	    print O "$boid\t$refseq\t$refscore\tNA\tNA\tNA\tNA\n";
	}
	else {
	    print O "$boid\t$refseq\t$refscore\t$ref->{$refseq}->{'start'}\t$ref->{$refseq}->{'end'}\t$ref->{$refseq}->{'strand'}\t$ref->{$refseq}->{'desc'}\n";
	    print ANNOT "$chrom\t$ref->{$refseq}->{'start'}\t$ref->{$refseq}->{'end'}\n";

	    # to check core outside hybrid region:
#	    if (($ref->{$refseq}->{'start'} < 2693400) and ($ref->{$refseq}->{'start'} > 382470)){
#		print "$boid\t$refseq\t$refscore\t$ref->{$refseq}->{'start'}\t$ref->{$refseq}->{'end'}\t$ref->{$refseq}->{'strand'}\t$ref->{$refseq}->{'desc'}\n";
#	    }
	    # to check core inside hybrid region:
	    if (($ref->{$refseq}->{'start'} > 2693400) or ($ref->{$refseq}->{'start'} < 382470)){
		print "$boid\t$refseq\t$refscore\t$ref->{$refseq}->{'start'}\t$ref->{$refseq}->{'end'}\t$ref->{$refseq}->{'strand'}\t$ref->{$refseq}->{'desc'}\n";
            }

	}
    }
    else {
	print O "$boid\tNONE\tNA\tNA\tNA\tNA\tNA\n";
    }
}
close (O);    
close (ANNOT);

###SUBS###

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
