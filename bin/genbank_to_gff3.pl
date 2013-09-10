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

my ($help, $gb);
GetOptions(
    'h|help'          => \$help,
    'b|gb=s'          => \$gb,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($gb){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

#####MAIN#####

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
    
}
