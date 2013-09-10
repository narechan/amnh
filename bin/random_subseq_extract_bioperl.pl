#!/usr/bin/perl -w

=head1 NAME

random_subseq_extract.pl

=head1 SYNOPSIS

random_subseq_extract.pl

Options:

--fasta is the fasta file that we need indexed
--number is the number of fasta segments we want
--length is the length of those segments
--outdir is your output dir

Requires the bioperl libs. 

=head1 DESCRIPTION

This program simulates reads from a source fasta file
and does so using both strands.

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

my ($help, $fasta, $number, $length, $outdir);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'n|number=s'      => \$number,
    'l|length=s'      => \$length,
    'o|outdir=s'      => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($fasta, $number, $length, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# read fasta file and count the number of seqs it holds                           
my $seqin = Bio::SeqIO -> new (-format => 'Fasta', -file => "$fasta");

my $counter = 0;
my %lengths;
my %ids;
while (my $sequence_obj = $seqin -> next_seq()){
    $counter++;
    my $id         = $sequence_obj->display_id();
    my $longness   = $sequence_obj->length();
    $lengths{$id}  = $longness;
    $ids{$counter} = $id;
}

# index the fasta file
my $index = Bio::Index::Fasta->new(-filename => $fasta . ".idx", -write_flag => 1);
$index->make_index($fasta);

# get the fasta file name
my $fasta_name;
if ($fasta =~/\//g){
    $fasta =~m/.*\/(.*)$/;
    $fasta_name = $1;
}

else {
    $fasta_name = $fasta;
}


# generating the random set
# NOTE: no checking for overlapping regions: OVERLAPS ALLOWED!  
open (INDX, ">$outdir/$fasta_name.randomindex");
open (FASTA, ">$outdir/$fasta_name.randomfasta");
my %random;
for (my $i = 1; $i <= $number; $i++){
    my $randchrom = int(rand($counter)) + 1;
    my $strand = int(rand(2)) + 1; # 1 is forward; 2 is reverse
    
    my $randstrand;
    if ($strand == 1){
	$randstrand = 1;
    }
    else{
	$randstrand = -1;
    }
    
    my $start;
    my $end = int rand($lengths{$ids{$randchrom}}) + 1;
    my $location;
    if ($randstrand == 1){
	$start = $end - $length + 1; 
	
	# check to see if the start is less than 1
	# POSSIBLE EDGE EFFECTS!
	if ($start < 1){
	    $i--;
	    next;
	}

	else{
	    $location = Bio::Location::Simple->new(-start  => $start,
						   -end    => $end,
						   -strand => $randstrand);
	}
	
    }
    else {
	$start = $end + $length - 1;
	
	# check to see if the start is greater than the length
	# POSSIBLE EDGE EFFECTS!
	if ($start > $lengths{$ids{$randchrom}}){
	    $i--;
	    next;
	}
	else{
	    $location = Bio::Location::Simple->new(-start  => $end,
						   -end    => $start,
						   -strand => $randstrand);
	}
	
    }
    
    my $sequence = $index->fetch($ids{$randchrom});
    my $subseq   = $sequence->subseq($location);

    if ($randstrand == 1){
        print FASTA ">$i:$ids{$randchrom}:$start:$end:$randstrand\n$subseq\n";
	print INDX "$randchrom\t$randstrand\t$start\t$end\n";
    }
    else {
        print FASTA ">$i:$ids{$randchrom}:$end:$start:$randstrand\n$subseq\n";
	print INDX "$randchrom\t$randstrand\t$end\t$start\n";
    }
    
}
close (INDX);
close (FASTA);
