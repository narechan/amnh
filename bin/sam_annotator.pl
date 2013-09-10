#!/usr/bin/perl -w

=head1 NAME

sam_annotator.pl

=head1 SYNOPSIS

sam_annotator.pl

Options:

--gff3 is your gff3 annotation file (optional)
--gb is your genbank annotation file (optional)
--sam is your sam file containing the alignments
--len is the length of your reads
--qthreshold is the mapping quality you're willing to accept
--verbose is set if you want read by read stats and final stats
    print to stderr

Requires the bioperl libs. 

=head1 DESCRIPTION

This program annotates the alns in a sam file given the annotations
in an annotation file. Either a gff3 or a genbank file must be supplied.

NOTE: GFF3/GB PROCESSING IS SPECIFIC TO BACTERIA:
    KEEPING CDS (MRNA), TRNA, RRNA, AND TRANSCRIPTS)

TODO: NEED TO UPDATE GFF PARSER TO PRIME THE DATA STRUC WITH GENE_IDS!!
    GB PARSER CURRENTLY DOES THIS
    
    CURRENTLY ALLOWS PARTIALLY OVERLAPPING READS AS COUNTS; SHOULD BE AN OPTION
    
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
use Bio::Seq;
use Bio::SeqIO;

my ($help, $gff3, $sam, $gb, $len, $qthresh, $verb);
GetOptions(
    'h|help'          => \$help,
    'g|gff3=s'        => \$gff3,
    's|sam=s'         => \$sam,
    'b|gb=s'          => \$gb,
    'l|len=s'        => \$len,
    'q|qthreshold=s'       => \$qthresh,
    'v|verbose'      => \$verb,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($sam, $len, $qthresh){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

# global data keepers
my $data = {};
my $cache = {};

# parse and store the annotations
my $gff3str = {};
my $gff3att = {};
my $gff3typ = {};

if ($gff3){ ## NEED TO UPDATE THIS TO PRIME THE DATA STRUC WITH GENE_IDS!!
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

        # assign a name and desc and pt id
	my $name;
	my $desc;
	my $protein_id;
        if ($feat->has_tag("gene")){
	    my @names = $feat->get_tag_values("gene");
            $name = join "_", @names;
        }
        if ($feat->has_tag("product")){
            my @descs = $feat->get_tag_values("product");
            $desc = join "_", @descs;
        }
	if ($feat->has_tag("protein_id")){
	    my @ids = $feat->get_tag_values("protein_id");
            $protein_id = join "_", @ids;
	}

	# get other crap                                                                                  
        my $start = $feat->start;
        my $end   = $feat->end;
        my $strand = $feat->strand;
        my $range = $start . "-" . $end;
	my $len = $end - $start + 1;
	
        if ($strand == 1){
            $gff3str->{$chrom}->{$range} = "+";
        }
        else {
            $gff3str->{$chrom}->{$range} = "-";
        }
	
	# cache the gene names and store as an attribute
	if ($name){
	    $cache->{$name} = $len;
	    $gff3att->{$chrom}->{$range} = $name;
	}
	elsif ($protein_id){
	    $cache->{$protein_id} = $len;
	    $gff3att->{$chrom}->{$range} = $protein_id;
	}
	else {
	    $cache->{$desc} = $len;
	    $gff3att->{$chrom}->{$range} = $desc;
        }

	$gff3typ->{$chrom}->{$range} = $type;
    }
}

else {
    print STDERR "Need an annotation file\n";
    die;
}

# cycle through the sam file
my $unmapped = 0;
my $lowmapped = 0;
my $mapped = 0;
open (S, "$sam");
while (my $line = <S>){
    next if ($line =~m/^\@/);
    chomp $line;

    my ($qryname,
	$flag,
        $refname,
        $refposst,
        $mapq,
        $cigar,
        $rnext,
        $pnext,
        $tlen,
        $seq,
        $qual,
	$alninfo) = split (/\t/, $line);

    # bailers: bail if no alignment or poor alignment
    if ($flag == 4){
	$unmapped++;
	next;
    }
    if ($mapq < $qthresh){
	$lowmapped++;
	next;
    }
    
    # otherwise look at the mapped read
    $mapped++;
    
    # find the end position of the read
    my $refposend = $refposst + $len;

    # find what kind of sequence the aln is in
    my $marker = 0;
    foreach my $range (sort keys %{$gff3str->{$refname}}){
	my ($start, $end) = split (/\-/, $range);
	
	# get attr
	my $gene = $gff3att->{$refname}->{$range};

	# get the strand                                                                                     
	my $strand;
	if ($gff3str->{$refname}->{$range} eq "+"){
	    $strand = 1;
	}
	else {
	    $strand = -1;
	}
	
	# search overlaps for the forward strand
	if ($strand == 1){
	    
	    # condition1: forward strand hanging over 5' end
	    if (($refposst < $start) and ($refposend >= $start)){
		$data->{$gene}->{'f5'}++;
		(print STDERR "$qryname\t$flag\t$refposst\t$refposend\t$gene\t$gene\t$end\tf5\n") if ($verb);
	    }
	    
	    # condition2: forward strand hanging over 3' end
	    elsif (($refposst <= $end) and ($refposend > $end)){
		$data->{$gene}->{'f3'}++;
		print STDERR ("$qryname\t$flag\t$refposst\t$refposend\t$gene\t$start\t$end\tf3\n") if ($verb);
	    }
	    
	    # condition3: inside the interval
	    elsif (($refposst >= $start) and ($refposend <= $end)){
		$data->{$gene}->{'fIN'}++;
		print STDERR ("$qryname\t$flag\t$refposst\t$refposend\t$gene\t$start\t$end\tfIN\n") if ($verb);
	    }
	    else {
		next;
	    }
	}
	
	# search overlaps for the reverse strand
	else {
	    
	    # condition1: reverse strand hanging over 3' end                                              
            if (($refposst < $start) and ($refposend >= $start)){
                $data->{$gene}->{'r3'}++;
                print STDERR ("$qryname\t$flag\t$refposst\t$refposend\t$gene\t$start\t$end\tr3\n") if ($verb);
            }

            # condition2: reverse strand hanging over 5' end                                         
            elsif (($refposst <= $end) and ($refposend > $end)){
                $data->{$gene}->{'r5'}++;
                print STDERR ("$qryname\t$flag\t$refposst\t$refposend\t$gene\t$start\t$end\tr5\n") if ($verb);
            }

            # condition3: inside the interval                                                             
            elsif (($refposst >= $start) and ($refposend <= $end)){
                $data->{$gene}->{'rIN'}++;
                print STDERR ("$qryname\t$flag\t$refposst\t$refposend\t$gene\t$start\t$end\tfIN\n") if ($verb);
            }
            else {
                next;
            }
	    
	}
    }
}

# calculate rpkm and output table
foreach my $gene (keys %$cache){
    if (exists ($data->{$gene})){
	my $counter = 0;
	foreach my $alntype (keys %{$data->{$gene}}){
	    $counter += $data->{$gene}->{$alntype};
	}
	my $rpkm = (($counter / ($cache->{$gene} / 1000)) / ($mapped / 1000000));
	print "$gene\t$counter\t$rpkm\n";
    }
    else {
	print "$gene\tNONE\n";
    }
}

if ($verb){
    print STDERR "MAPPED: $mapped\n";
    print STDERR "LOWQUAL: $lowmapped\n";
    print STDERR "UNMAPPED: $unmapped\n";
}
