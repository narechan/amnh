#!/usr/bin/perl

use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

my %opts = ();
getopts ('f:o:n:', \%opts);
my $fasta  = $opts{'f'};
my $outdir = $opts{'o'};
my $number = $opts{'n'};

`mkdir -p $outdir`;

# read fasta file and count the number of seqs it holds
my $seqin = Bio::SeqIO -> new (-format => 'Fasta', -file => "$fasta");

my $counter = 0;
print stderr "Counting up the seqs\n";
while (my $sequence_obj = $seqin -> next_seq()){
    $counter++;
}

# find the number to put in each fasta file
my $switcher = int (($counter / $number) + 0.5);

# build the subfastas
my $iterator = 0;
my $multiple = 1;

print stderr "Writing subfasta $multiple\n";
open (fh, ">$outdir/query.$multiple");
my $sequence_in = Bio::SeqIO -> new (-format => 'Fasta', -file => "$fasta");
while (my $sequence_obj = $sequence_in -> next_seq()){
    $iterator++;  
    my $id       = $sequence_obj -> display_id();
    my $sequence = $sequence_obj -> seq();
    my $desc     = $sequence_obj->desc();
    
    if ($iterator <= ($switcher * $multiple)){
	print fh ">$id $desc\n$sequence\n";
    }

    # add the remainder to the end of the last file 
    elsif ((($switcher * $number) < $iterator) and ($iterator <= $counter)){
	print fh ">$id $desc\n$sequence\n";
    }
    else {
	close (fh);
	$multiple++;

	print stderr "Writing subfasta $multiple\n";
	
	open (fh, ">$outdir/query.$multiple");
	print fh ">$id $desc\n$sequence\n";
    }
}

