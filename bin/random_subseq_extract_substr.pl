#!/usr/bin/perl -w

=head1 NAME

random_subseq_extract.pl

=head1 SYNOPSIS

random_subseq_extract.pl

Options:

--fasta is the fasta file
--number is the number of fasta segments we want
--length is the length of those segments
--outdir is your output dir

Requires the bioperl libs. 

NOTE THAT THE SEQUENCES MUST BE ON THE SAME LINE IN THE READS FILE

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
use Bio::Seq;
use Bio::SeqIO;
use Bio::Index::Fasta;
use Bio::Location::Simple;

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
open (F, "$fasta");
my $counter = 0;
my %lengths;
my %ids;
my %seqs;
my %revseqs;

my $id;
my $seq;
my $len;
my $revseq;
while (my $line = <F>){
    chomp $line;

    if ($line =~m/^>(.*)/){
        $id = $1;
        $counter++;
	$ids{$counter} = $id;
    }
    else {
        $seq = $line;
	$revseq = revdnacomp ($seq);
	$seqs{$id} = $seq;
        $revseqs{$id} = $revseq;

        $len = length ($seq);
        $lengths{$id} = $len;
    }
}
close (F);

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
    my $subseq;
    if ($randstrand == 1){
	$start = $end - $length + 1; 
	
	# check to see if the start is less than 1
	# POSSIBLE EDGE EFFECTS!
	if ($start < 1){
	    $i--;
	    next;
	}
	else{
	    $subseq = substr($seqs{$ids{$randchrom}}, $start - 1, $length); 
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
	    $subseq = revdnacomp (substr($seqs{$ids{$randchrom}}, $end - 1, $length));
	}
    }
    
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

####SUBS####

sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);

    $revcomp =~ tr/ACGTacgt/TGCAtgca/;

    return $revcomp;
}
