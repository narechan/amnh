#!/usr/bin/perl

### parse_fasta1.pl

## This program accepts a fasta file as an arguement,
## parses the fasta file, and computes the min, max and 
## average lengths of the sequences it contains.

# take first element off the perl special array variable 
my $fastafile = shift @ARGV;

# declare global variables to track total sequence length,
# count the number of sequences, and track max and min lengths
my $totalSeqLen = 0;
my $seqCounter  = 0;
my $max         = 0;
my $min         = 1e100;

# open the fastafile and go through it line by line
open (FASTA, "$fastafile");            
while (my $line = <FASTA>){            
    chomp $line;                       
    
    # declare local variables to track sequence id, sequence
    # and sequence length
    my $id;
    my $sequence;
    my $seqLen;
    
    # parse the defline
    if ($line =~m/^>/){
	$seqCounter++;
	$id = $line;
    }

    # parse the sequence
    else {
	$sequence = $line;
	
	# get sequence length and update total
	# sequence length
	$seqLen = length ($sequence);
	$totalSeqLen += $seqLen;
	
	# update min and max
	if ($seqLen > $max){
	    $max   = $seqLen;
	}
	if ($seqLen < $min){
	    $min   = $seqLen;
	}
    }
}

# compute the average
my $avgLen = $totalSeqLen / $seqCounter;

# print your report
print "Avg = $avgLen\n";
print "Max = $max\n";
print "Min = $min\n";
