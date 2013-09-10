#!/usr/bin/perl -w

=head1 NAME

insidP.pl

=head1 SYNOPSIS

The script generates files of reads containing incrementally 
increasing random fractions of sequence from a total read set given the 
increment (e.g. tenths) or a file giving the desired reads per subset,
or a file with the total number of bases per subset.

Note that the files can contain repeated reads (bootstrap -- sample with
replacement) or can contain no repeated reads (jackknife -- sample without 
replacement).

Note that either factor, reads, or bases must be chosen to quantify 
reads desired.

NOTE THAT THE SEQUENCES MUST BE ON A SINGLE LINE IN THE INPUT READS FILE!

Dependencies:                                                                                            

None.

Options:

--infile is the input reads file
--outdir is the output directory
--mode is the sampling mode (bootstrap = 1; jackknife = 2)
--factor is the sequence fraction you 
    want to increment up to the 
    entire dataset (optional)
--reads is a file containg the read counts
    per subset (optional)
--bases is a file containing the base counts
    per subset (optional)

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2010 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $infile, $outdir, $factor, $readfile, $basefile, $mode);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'o|outdir=s'      => \$outdir,
    'f|factor=s'      => \$factor,
    'r|reads=s'       => \$readfile,
    'b|bases=s'       => \$basefile,
    'm|mode=s'        => \$mode,
    ) or pod2usage;
pod2usage if $help;

for my $option ($infile, $outdir, $mode){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####    

# read fasta file and count the number of seqs it holds                                                 
open (F, "$infile");
my $counter = 0;
my $totlen  = 0;
my $hashlen = {};
while (my $line = <F>){
    chomp $line;

    my $id;
    my $seq;
    if ($line =~m/^>/){
	my @line = split (/\s/, $line);
	$id = $line[0];
	$counter++;
    }
    else {
	$seq = $line;
	my $len = length ($seq);
        $totlen += $len;
        $hashlen->{$counter} = $len;
    }
}
close (F);

# find the number to put in each fasta file                                                       
my @subfastas;
my @subcounts;
if ($factor){
    my $frac = int (($counter * $factor) - 0.5);
    my $seqcount = $frac;
    until ($seqcount >= $counter){
	push (@subfastas, $seqcount);
	$seqcount += $frac;
    }
    # pop off last value and put in total to acct for remainders                                          
    pop @subfastas;
    push (@subfastas, $counter);
}
elsif ($readfile){
    open (RF, "$readfile");
    while (my $line = <RF>){
	chomp $line;
	
	if ($line > $counter){
	    print STDERR "Not enough seqs for your subset\n";
	    die;
	}
	else {
	    push (@subfastas, $line);
	}
    }
    close (RF);
}
elsif ($basefile){
    open (BF, "$basefile");
    while (my $line = <BF>){
	chomp $line;
	
	if ($line > $totlen){
	    print STDERR "Not enough seqs for your basecount\n";
	    die;
	}
	else{
	    push (@subcounts, $line);
	}
    }
    close (BF);
}
else {
    print STDERR "reads, bases, or factor must be defined\n";
    die;
}

# build the subfastas 
if ($readfile or $factor){
    foreach my $subfasta (@subfastas){
	print STDERR "Creating $subfasta\n";
	open (OUT, ">$outdir/$subfasta.fa");
	
	my %population;
	for (my $sample_el = 1; $sample_el <= $subfasta; $sample_el++){
	    my $random_idx = int(rand($counter)) + 1;
	    if ($mode == 2){
		if (exists ($population{$random_idx})){
		    ($random_idx = int(rand($counter)) + 1) until
			(! exists($population{$random_idx}));
		}
	    }
	    $population{$random_idx}++;
	}
	
	my $subcounter = 0;
	my $id;
	my $seq;
	open (FH, "$infile");
	while (my $line = <FH>){
	    chomp $line;
	    
	    if ($line =~m/^>/){
		my @line = split (/\s/,$line);
		$id = $line[0];
	    }
	    else {
		$seq = $line;
		$subcounter++;
		
		if (exists ($population{$subcounter})){
		    for (my $sample = 1; $sample <= $population{$subcounter}; $sample++){
                        print OUT "$id.$sample\n$seq\n";
                    }
		}
		else {
		    next;
		}
	    }
	}
	close (OUT);
    }
}
else {
    foreach my $subcount (@subcounts){
	print STDERR "Creating $subcount\n";
	open (OUT, ">$outdir/$subcount.fa");
	
	my %population;
	my $length = 0;
	until ($length > $subcount){
	    my $random_idx = int(rand($counter)) + 1;
	    if ($mode == 2){
		if (exists ($population{$random_idx})){
		    ($random_idx = int(rand($counter)) + 1) until
			(! exists($population{$random_idx}));
		}
	    }
	    $population{$random_idx}++;
	    $length += $hashlen->{$random_idx};
	}
	    
	my $subcounter = 0;
	my $id;
	my $seq;
	open (FH, "$infile");
	while (my $line = <FH>){
	    chomp $line;
	    
	    if ($line =~m/^>/){
		my @line = split (/\s/,$line);
		$id = $line[0];
	    }
	    else {
		$seq = $line;
		$subcounter++;
		
		if (exists ($population{$subcounter})){
		    for (my $sample = 1; $sample <= $population{$subcounter}; $sample++){
			print OUT "$id.$sample\n$seq\n";
		    }
		}
		else {
		    next;
		}
	    }
	}
	close (OUT);
    }
}
