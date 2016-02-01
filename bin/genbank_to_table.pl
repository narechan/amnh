#!/usr/bin/perl -w

=head1 NAME

genbank_to_gff3.pl.pl

=head1 SYNOPSIS

Options:

--gb is your genbank annotation file

Requires the bioperl libs. 

=head1 DESCRIPTION

GB --> GFF3
    
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
use Bio::Index::Fasta;

my ($help, $gb, $outdir, $reffa);
GetOptions(
    'h|help'          => \$help,
    'b|gb=s'          => \$gb,
    'o|outdir=s'      => \$outdir,
    'r|reffa=s'       => \$reffa,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($gb, $outdir, $reffa){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

print STDERR "Processing the reference\n";
# index the fasta file                                                              
my $index = Bio::Index::Fasta->new(-filename => $reffa . ".idx", -write_flag => 1);
$index->make_index($reffa);

my $seqin = Bio::SeqIO->new(-file=>$gb, -format=>'Genbank');       
my $seqobj = $seqin->next_seq;
my $chrom  = $seqobj->display_id;

# cycle through the features
open (P, ">$outdir/pt.fa");
open (T, ">$outdir/table");
foreach my $feat ($seqobj->get_SeqFeatures){
    my $type = $feat->primary_tag;
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
    
    print T "$name\t$start\t$end\t$strand\t$desc\n";
}
close (P);
close (T);
