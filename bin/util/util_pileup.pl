#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($help, $infile, $outdir);
GetOptions(
    'h|help'           => \$help,
    'i|infile=s'       => \$infile,
    'o|outdir=s'       => \$outdir,	   
    ) or pod2usage;

pod2usage if $help;

# get file name                                                                                     
my $infile_name;
if ($infile =~/\//g){
    $infile =~m/.*\/(.*)$/;
    $infile_name = $1;
}

else {
    $infile_name = $infile;
}


# parse the pileup
open (PU, "$infile");
open (MI, ">$outdir/$infile_name.snpsonly");
open (QF, ">$outdir/$infile_name.snpsfiltered");
open (ID, ">$outdir/$infile_name.indels");
while (my $line = <PU>){
    chomp $line;
    my ($refseq, $refpos, $refbase, $snp, $cq, $sq, $rms, $cov, $bases, $quals, $rest) 
	= split (/\t/, $line);
    if ($refbase eq "*"){
	print ID "$line\n";
    }
    else{
	if ($sq > 0){
	    print MI "$line\n";
	}
	else {
	    next;
	}
	
	if ($sq >= 20){
	    print QF "$line\n";
	}
	else {
	    next;
	}
    }
}
close (ID);
close (QF);
close (MI);
close (PU);
