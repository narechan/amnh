#!/usr/bin/perl

#####                                                                             
#####                                                                             

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $infile, $refgb);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'g|refgb=s'       => \$refgb,
	   ) or pod2usage;

pod2usage if $help;

# parse the genbank file
my $ref = {};
my $seqin = Bio::SeqIO->new(-file=>$refgb, -format=>'Genbank');
my $seqobj = $seqin->next_seq;
my $chrom  = $seqobj->display_id;

# cycle through the features                                                 
foreach my $feat ($seqobj->get_SeqFeatures){
    my $type = $feat->primary_tag;
    
    next unless ($feat->primary_tag =~m/CDS|tRNA|rRNA|transcript/);

    # parse the data
    my $gene;
    my $name;
    my $prod;
    if ($feat->has_tag("gene")){
	my @names = $feat->get_tag_values("gene");
	$gene = join "_", @names;
    }
    else {
	$gene = "NA";
    }
    if ($feat->has_tag("locus_tag")){
	my @names = $feat->get_tag_values("locus_tag");
	$name = join "_", @names;
    }
    else {
	print STDERR "NO TAG!\n";
	die;
    }
    if ($feat->has_tag("product")){
	my @names = $feat->get_tag_values("product");
	$prod = join "_", @names;
    }
    else {
	$prod = "NA";
    }
    
    my $start = $feat->start;
    my $end   = $feat->end;
    my $strand = $feat->strand;
    
    # store data
    $ref->{$name}->{'start'} = $start;
    $ref->{$name}->{'end'} = $end;
    $ref->{$name}->{'strand'} = $strand;
    $ref->{$name}->{'prod'} = $prod;
    $ref->{$name}->{'gene'} = $gene;
}

open (I, "$infile");
while (my $line = <I>){
    next if ($line =~m/^\"\"/);
    chomp $line;
    my @fields = split (/\,/, $line);
    my $acc = $fields[1];
    $acc =~s/\"//g;
    
    print "$line,$ref->{$acc}->{'gene'},$ref->{$acc}->{'prod'}\n";
}
close (I);

