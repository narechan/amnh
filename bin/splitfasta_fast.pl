#!/usr/bin/perl -w


#####SETUP#####
use strict;
use Pod::Usage;
use Getopt::Long;

my ($help, $fasta, $number, $outdir);
GetOptions(
    'h|help'          => \$help,
    'f|fasta=s'       => \$fasta,
    'n|number=s'      => \$number, #number of seqs in each file
    'o|outdir=s'      => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($fasta, $number, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# count the number of seqs
print STDERR "Counting up the seqs in your fasta file\n";
open (A, "$fasta");
my $seqcounter = grep /^>/, <A>;
close (A);

#my $seqcounter = 0;
#open (A, "$fasta");
#while (my $line = <A>){
#    chomp $line;
#    print STDERR "$seqcounter\n";
#    if ($line =~m/^>(.*)/){
#        $seqcounter++;
#    }
#    else {
#	next;
#    }
#}
#close (A);

# find the number to put in each fasta file                                          
my $switcher = int (($seqcounter / $number) + 0.5);

# iterate
open (F, "$fasta");
my $counter = 0;
my $totcounter = 0;
my $multiple = 1;
my $id;
my $seq;
open (FH, ">$outdir/query.$multiple");
open (F, "$fasta");
while (my $line = <F>){
    chomp $line;

    if ($line =~m/^>(.*)/){
        $id = $1;
        $counter++;
	$totcounter++;
    }
    else {
        $seq = $line;
	
	if ($totcounter == $seqcounter){
	    print FH ">$id\n$seq\n";
	    close (FH);
	}
	elsif ($counter > $switcher){
	    close (FH);
	    $multiple++;
	    
	    print STDERR "Writing subfasta $multiple\n";
	    
	    open (FH, ">$outdir/query.$multiple");
	    print FH ">$id\n$seq\n";
	    $counter = 1;
	}
	else {
	    print FH ">$id\n$seq\n";
	}
    }
}
close (F);
