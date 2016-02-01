#!/usr/bin/perl -w

=head1 NAME

snp_annotator.pl

=head1 SYNOPSIS

snp_annotator.pl

Options:

--gff3 is your gff3 annotation file (optional)
--gb is your genbank file (optional)
--fasta is your contig / chrom fasta file
--vcf is your vcf file containing variant information (optional)
--table is a nucmer table of snps for finished seqs (optional)
--config is your config file
    (contains the codon table of your choice)
--quality is the quality score cutoff for reporting (optional)
--query is your query name (for labeling purposes)    

Requires the bioperl libs. 

=head1 DESCRIPTION

This program annotates snps in a vcf file given the annotations
in an annotation file. Either a ggf3 or a genbank file must be supplied.
Either a vcf file or nucmer table must be supplied

NOTE: ONLY PROCESSES SNPS IN VCF FOR NOW (NO INDELS)
NOTE: GFF3 PROCESSING IS SPECIFIC TO BACTERIA:
    KEEPING CDS (MRNA), TRNA, RRNA, AND TRANSCRIPTS)
NOTE: SKIPPING CODON ANALYSIS OF HETEROZYGOTES (spurious for bacteria)

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

#####SETUP#####
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Index::Fasta;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $fasta, $gff3, $vcf, $config, $query, $gb, $table, $quality);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'g|gff3=s'        => \$gff3,
    'v|vcf=s'         => \$vcf,
    'c|config=s'      => \$config,	   
    'q|query=s'       => \$query,
    'b|gb=s'          => \$gb,
    't|table=s'       => \$table,	   
    'x|quality=s'     => \$quality,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($fasta, $config, $query){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# process the config file
my $configuration = load_config ($config);

# index the fasta file
my $index = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$index->make_index($fasta);

# store the codon table
my $codons = {};
my $codontable = $configuration->{'CODON_TABLE'};
my @codons = split (/\,\s*/, $codontable);
foreach my $codon (@codons){
    my ($cd, $aa) = split (/\//, $codon);
    $codons->{$cd} = $aa;
}

# parse and store the annotations
my $gff3str = {};
my $gff3att = {};
my $gff3typ = {};

if ($gff3){
    open (G, "$gff3");
    while (my $line = <G>){
	next if ($line =~m/^\#/);
	chomp $line;
	
	my ($chrom,
	    $source,
	    $type,
	    $start,
	    $end,
	    $score,
	    $strand,
	    $phase,
	    $attr) = split (/\t/, $line);
	
	my @attrs = split (/\;/, $attr);
	
	# store the information for tRNA, rRNA, CDS (mRNA), 
	# and transcript only
	if ($type =~m/CDS|tRNA|rRNA|transcript/){
	    my $range = $start . "-" . $end;
	    $gff3str->{$chrom}->{$range} = $strand;
	    $gff3att->{$chrom}->{$range} = [@attrs];
	    $gff3typ->{$chrom}->{$range} = $type;
	}
	else {
	    next;
	}
    }
    close (G);
}

elsif ($gb){
    my $seqin = Bio::SeqIO->new(-file=>$gb, -format=>'Genbank');       
    my $seqobj = $seqin->next_seq;
    my $chrom  = $seqobj->display_id;

    # cycle through the features                                                                      
    foreach my $feat ($seqobj->get_SeqFeatures){
	my $type = $feat->primary_tag;
        next unless ($feat->primary_tag =~m/CDS|tRNA|rRNA|transcript/);

        # assign a name and desc
        my @attrs;
	my $name;
	my $desc;
        if ($feat->has_tag("gene")){
	    my @names = $feat->get_tag_values("gene");
            $name = join "_", @names;
	    push (@attrs, "Name=$name");
        }
        if ($feat->has_tag("product")){
            my @descs = $feat->get_tag_values("product");
            $desc = join "_", @descs;
	    push (@attrs, "description=$desc");
        }
	my $start = $feat->start;
	my $end   = $feat->end;
	my $strand = $feat->strand;
	my $range = $start . "-" . $end;
	
	if ($strand == 1){
	    $gff3str->{$chrom}->{$range} = "+";
	}
	else {
	    $gff3str->{$chrom}->{$range} = "-";
	}

	$gff3att->{$chrom}->{$range} = [@attrs];
	$gff3typ->{$chrom}->{$range} = $type;
    }
}

else {
    print STDERR "Need an annotation file\n";
    die;
}

# print the header
print "Chrom\tQuery\tRefPos\tRefBase\tSNP\tSnpQual\tState\tRefCodon\tSNPCodon\tType\tName\tDesc\n";

# cycle through the vcf file and annotate only the snps    
if ($vcf){
    open (V, "$vcf");
}
else {
    open (V, "$table");
}

while (my $line = <V>){
    next if ($line =~m/^\#/);
    chomp $line;

    my ($refseq,
	$qry,
        $refpos,
        $rsid,
        $refbase,
        $snp,
        $sq,
        $filter,
        $info,
        $format,
        $genotypeinfo);

    if ($vcf){
	($refseq, 
	 $refpos, 
	 $rsid, 
	 $refbase, 
	 $snp, 
	 $sq, 
	 $filter, 
	 $info, 
	 $format, 
	 $genotypeinfo) = split (/\t/, $line);
    }
    else {
	($refseq,
	 $qry,
	 $refpos,
	 $refbase,
	 $snp,
	 $sq) = split (/\t/, $line);
    }

    # if quality cutoff is defined, apply here
    if ($quality){
	next if ($sq < $quality);
    }
    
    # skip all indels in vcf files
    if ($vcf){
	next if ($info =~m/INDEL/);
	next if ($snp =~m/^I/); #this deals with pacbio vcf files which are version 3.3
	next if ($snp =~m/^D/); #same
    }
    
    # find what kind of sequence the variant is in
    my $marker = 0;
    foreach my $range (sort keys %{$gff3str->{$refseq}}){
	my ($start, $end) = split (/\-/, $range);
	
	# if it's in a gene, find out which one, pull the sequence,
	# and annotate as syn or nonsyn
	if (($start <= $refpos) and ($refpos <= $end)){

	    # check to see whether vcf indicates two alleles;
	    # if so, report and skip
	    if ($snp =~m/\,/){
		print "$refseq\t$query\t$refpos\t$refbase\t$snp\t$sq\theterozygous\n";
		$marker = 1;
		last;
	    }

	    # check to see if the reference is ambiguous;
	    # if so, report and skip
	    if ($refbase eq "N"){
		print "$refseq\t$query\t$refpos\t$refbase\t$snp\t$sq\tambiguous\n";
		$marker = 1;
		last;
	    }

	    # get the gene annotations -- selected vars are empty if 
	    # key/val pairs don't exist in the attr
	    my @attrs = @{$gff3att->{$refseq}->{$range}};
	    my $name;
	    my $desc;
	    foreach my $attr (@attrs){
		my ($key, $value) = split (/\=/, $attr);
		if ($key eq "Name"){
		    $name = $value;
		}
		if ($key eq "description"){
		    $desc = $value;
		}
	    }
	    
	    # get the strand
	    my $strand;
	    if ($gff3str->{$refseq}->{$range} eq "+"){
		$strand = 1;
	    }
	    else {
		$strand = -1;
	    }

	    # get the sequence
	    my $location = Bio::Location::Simple->new(-start  => $start,
						      -end    => $end,
						      -strand => $strand);
	    my $sequence  = $index->fetch($refseq);
	    my $subseq    = $sequence->subseq($location);
	    my $revsubseq = reverse ($subseq);
	    
	    # codon specific tracker variables
	    my $state;
            my $refstate;
            my $altstate;
            my $st = $start;
            my $ed;

	    # operate on the forward strand in a standard way
	    if ($strand == 1){
		my @codons   = $subseq =~m/(...)/g;
		
		foreach my $codon (@codons){
		    my $ed = $st + 2;
		    if (($st <= $refpos) and ($refpos <= $ed)){
			my $refaa = $codons->{$codon};
			$refstate = $codon . "/" . $refaa;
			
			my $altcodon;
			my $altaa;
			my @codon = split (//, $codon);
			if ($st == $refpos){
			    $altcodon = $snp . $codon[1] . $codon[2];
			    $altaa    = $codons->{$altcodon};
			    $altstate = $altcodon . "/" . $altaa;
			}
			elsif ($ed == $refpos){
			    $altcodon = $codon[0] . $codon[1] . $snp;
			    $altaa    = $codons->{$altcodon};
			    $altstate = $altcodon . "/" . $altaa;
			}
			else {
			    $altcodon = $codon[0] . $snp . $codon[2];
			    $altaa    = $codons->{$altcodon};
			    $altstate = $altcodon . "/" . $altaa;
			}
			
			# to allow for stuff like 'V+' in the codon table
			if (($refaa =~m/\Q$altaa/) or ($altaa =~m/\Q$refaa/)){
			    $state = "synonymous";
			}
			else {
			    $state = "nonsynonymous";
			}
			last;
			
		    }
		    else {
			$st = $ed + 1;
			next;
		    }
		}
	    }

	    # for genes on the reverse strand with snps described with respect to 
	    # the forward strand, take the following steps: Operate on the flipped 
	    # reverse complement to put the gene string in the snp's forward coordinate
	    # system. Reverse the codons as analyzed. Reverse complement the snps as 
	    # given on the forward strand. Reverse the $st and $ed coordinates.
	    else {
		my @codons   = $revsubseq =~m/(...)/g;
		foreach my $codon (@codons){
		    my $revcodon = reverse ($codon);
		    my $revcompsnp = revdnacomp ($snp);
                    my $ed = $st + 2;
		    
                    if (($st <= $refpos) and ($refpos <= $ed)){
                        my $refaa = $codons->{$revcodon};
                        $refstate = $revcodon . "/" . $refaa;

                        my $altcodon;
                        my $altaa;
                        my @codon = split (//, $revcodon);
                        if ($ed == $refpos){
                            $altcodon = $revcompsnp . $codon[1] . $codon[2];
                            $altaa    = $codons->{$altcodon};
                            $altstate = $altcodon . "/" . $altaa;
                        }
                        elsif ($st == $refpos){
                            $altcodon = $codon[0] . $codon[1] . $revcompsnp;
                            $altaa    = $codons->{$altcodon};
                            $altstate = $altcodon . "/" . $altaa;
                        }
                        else {
                            $altcodon = $codon[0] . $revcompsnp . $codon[2];
                            $altaa    = $codons->{$altcodon};
                            $altstate = $altcodon . "/" . $altaa;
                        }

                        if (($refaa =~m/\Q$altaa/) or ($altaa =~m/\Q$refaa/)){
                            $state = "synonymous";
                        }
                        else {
                            $state = "nonsynonymous";
                        }
                        last;

                    }
                    else {
                        $st = $ed + 1;
                        next;
                    }
                }
		
	    }
	    if ($state){
		print "$refseq\t$query\t$refpos\t$refbase\t$snp\t$sq\t$state\t";
		print "$refstate\t$altstate\t$gff3typ->{$refseq}->{$range}\t";
	    }
	    else {
		print "$refseq\t$query\t$refpos\t$refbase\t$snp\t$sq\tNO STATE\t";
                print "NO STATE\tNO STATE\t$gff3typ->{$refseq}->{$range}\t";
	    }
	    if ($name){
		print "$name\t";
	    }
	    else{
		print "NO NAME\t";
	    }
	    if ($desc){
		print "$desc\n";
	    }
	    else {
		print "NO DESC\n";
	    }

	    $marker = 1;
#	    last;
	}
	else {
	    next;
	}
    }	
    
    # will designate as intergenic anything that did not satisfy the gene ranges;
    # even heterozygotes
    (print "$refseq\t$query\t$refpos\t$refbase\t$snp\t$sq\tintergenic\n") if 
	($marker == 0);
}
			  
	    
##### SUBS #####

sub load_config{
    my $file = shift;

    # get and store key/val pairs in config
    my $config = {};
    open (F, "$file") or
        die "Unable to open config file\n";
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) = 
            ($line =~ m/^\s*(\w+)\s*=\s*(.*)$/);
        $config->{$key} = $value;
    }
    close (F);

    return ($config);
}

sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);

    $revcomp =~ tr/ACGTacgt/TGCAtgca/;

    return $revcomp;
}
