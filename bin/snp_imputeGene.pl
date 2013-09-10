#!/usr/bin/perl -w

# This program accepts a list of snps from snp_annotator.pl and a 
# reference gff3 file and creates a set of imputed genes using the reference.

# NOTE: the number of taxa per alignment is the number of queries in 
#       your list plus 1, for the reference. so query list should not contain the reference

#       Also, this process only respects polymorphisms -- no data for indels

# -i is the input snp table
# -l is the query list
# -o is the outdir
# -q is the snp quality above which snps are accepted
# -g is the gff3 file (optional)
# -b is the gb file (optional)
# -s is the reference genome sequence

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Index::Fasta;
use Bio::Location::Simple;

my ($infile, $list, $outdir, $quality, $gff3file, $gbfile, $sequencefile);
GetOptions(
	   'i|infile=s'  => \$infile,
	   'l|list=s'    => \$list,
	   'o|outdir=s'  => \$outdir,
	   'q|quality=s' => \$quality,
	   'g|gff3file=s' => \$gff3file,
	   'b|gbfile=s'   => \$gbfile,
	   's|sequencefile=s' => \$sequencefile,
	   );

`mkdir -p $outdir/genes`;

#####MAIN#####

# read in the refenece sequence
my $refgenome = {};
my $refid;
my $seq;
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$sequencefile");
while (my $sequence_obj = $seqin->next_seq()){
    $refid  = $sequence_obj->id();
    $seq = $sequence_obj->seq();
}
my @seq = split (//, $seq);
my $poscounter = 0;
foreach my $seqpos (@seq){
    $poscounter++;
    $refgenome->{$poscounter} = $seqpos;
}

# read in the quey list
my @queries;
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    push (@queries, $line);
}
close (L);

# store the gff3 info
my $gff3str = {};
my $gff3att = {};
my $gff3typ = {};

if ($gff3file){
    open (G, "$gff3file");
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
elsif ($gbfile){
    my $seqin = Bio::SeqIO->new(-file=>$gbfile, -format=>'Genbank');
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
	elsif ($feat->has_tag("locus_tag")){
	    my @names = $feat->get_tag_values("locus_tag");
            $name = join "_", @names;
            push (@attrs, "Name=$name");
	}
	else {
	    print STDERR "no name!\n";
	    die;
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


# read in the snp file and store data
my $snps = {};
my $refbases = {};
open (I, "$infile");
while (my $line = <I>){
    chomp $line;
    my ($ref, $query, $refpos, $refbase, $snp, $sq, $stuff) = 
	split (/\t/, $line, 7);

    # plug in the reference chars
    $snps->{$refpos}->{$ref}   = $refbase;
    $refbases->{$refpos} = $refbase;
    
    # filters
    next if ($refbase eq "N"); #ambiguous ref
    if ($sq < $quality){ #sq
	$snps->{$refpos}->{$query} = "?";
	next;
    }
    if ($snp =~m/\,/){ #hets
	$snps->{$refpos}->{$query} = "?";
	next;
    }

    # plug in the snp chars
    $snps->{$refpos}->{$query} = $snp;
}
close (I);

# post process datastruc to fill in ref alleles
# for queries with snp gaps
foreach my $pos (keys %$snps){
    foreach my $query (@queries){
	if (exists($snps->{$pos}->{$query})){
	    next;
	}
	else {
	    $snps->{$pos}->{$query} = $refbases->{$pos};
	}
    }
}

# write out whole genome imputed polymorphism alignment                                                      
my $multgenaln = {};
foreach my $pos (sort {$a <=> $b} keys %$refgenome){
#    print STDERR "$pos\n";
    foreach my $query (@queries){
	if (exists($snps->{$pos}->{$query})){
	    $multgenaln->{$query} .= $snps->{$pos}->{$query};
	}
	else {
	    $multgenaln->{$query} .= $refgenome->{$pos};
	}
    }
    $multgenaln->{$refid} .= $refgenome->{$pos};
}

# print out the multgenerefid
open (GA, ">$outdir/wga.fa");
foreach my $taxon (sort keys %$multgenaln){
    print GA ">$taxon\n$multgenaln->{$taxon}\n";
}
close (GA);

# index the wga and take print out genes by looping though gff3 features
my $index = Bio::Index::Fasta->new(-filename => "$outdir/wga.fa.idx", -write_flag => 1);
$index->make_index("$outdir/wga.fa");

foreach my $range (sort keys %{$gff3str->{$refid}}){
    my ($start, $end) = split (/\-/, $range);
    
    # name
    my @attrs = @{$gff3att->{$refid}->{$range}};
    my $name;
    foreach my $attr (@attrs){
	my ($key, $value) = split (/\=/, $attr);
	if ($key eq "Name"){
	    $name = $value;
	}
    }
    print STDERR "$name\n";
    
    # strand
    my $strand;
    if ($gff3str->{$refid}->{$range} eq "+"){
	$strand = 1;
    }
    else {
	$strand = -1;
    }

    # sequence
    my $location = Bio::Location::Simple->new(-start  => $start,
					      -end    => $end,
					      -strand => $strand);
    
    open (OUT, ">$outdir/genes/$name.fa");
    foreach my $query (@queries){
	my $sequence  = $index->fetch($query);
	my $subseq    = $sequence->subseq($location);
	
	print OUT ">$query\n$subseq\n";
    }
    close (OUT);
}
