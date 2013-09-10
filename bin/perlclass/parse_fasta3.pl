#!/usr/bin/perl

### parse_fasta3.pl

## This program accepts a fasta file as an arguement,
## parses the fasta file, and computes the min, max and 
## average lengths of the sequences it contains.
## It also bins the lengths into a frequency distribution.

# take first element off the perl special array variable
my $fastafile = shift @ARGV;

# declare global variables to track total sequence length,
# count the number of sequences, track max and min lengths,
# and to compute bins
my $totalSeqLen = 0;
my $seqCounter  = 0;
my $max         = 0;
my $min         = 1e100;
my %bins;

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
	($max = $seqLen) if ($seqLen > $max);
	($min = $seqLen) if ($seqLen < $min);
	
	# update the bins
	if ($seqLen <= 400){
	    $bins{1}++;
	}
	elsif (($seqLen >= 401) and ($seqLen <= 800)){
	    $bins{2}++;
	}
	elsif (($seqLen >= 801) and ($seqLen <= 1200)){
	    $bins{3}++;
	}
	else {
	    $bins{4}++;
	}
    }
}

# compute the average
my $avgLen = $totalSeqLen / $seqCounter;

# print your report
print "Avg = $avgLen\n";
print "Max = $max\n";
print "Min = $min\n";
    
foreach my $bin (sort {$a <=> $b} keys %bins){
    print "Bin $bin = $bins{$bin}\n";
}
