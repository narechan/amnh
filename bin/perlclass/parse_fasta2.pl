#!/usr/bin/perl

### parse_fasta2.pl

## This program accepts a fasta file as an arguement,                                                     
## parses the fasta file, and computes the min, max and                                                    
## average lengths of the sequences it contains
## BUT DOES IT ALL USING EXTERNAL PERL LIBRARIES. 

# invoke the modules that you will be using
use Bio::Seq;
use Bio::SeqIO;
use Statistics::Descriptive;

# take first element off the perl special array variable
my $fastafile = shift @ARGV;

# declare a global array to store each sequence's length
my @lens;

# open the fasta file using bioperl's SeqIO module
# and iterate record by record
my $seqin = Bio::SeqIO->new (-format=>'Fasta', -file=>"$fastafile");
while (my $sequence_obj = $seqin->next_seq()){

    # get id, sequence, and sequence length
    my $id       = $sequence_obj->display_id();
    my $seq      = $sequence_obj->seq();
    my $len      = length ($seq);
	
    # store sequence lengths into the array
    push (@lens, $len);
}

# generate statistics on sequence lengths using the 
# Statistics::Descriptive module
my $statobj = Statistics::Descriptive::Full->new();
$statobj->add_data(@lens);
my $avg = $statobj->mean();
my $max = $statobj->max();
my $min = $statobj->min();

# print your report
print "Avg = $avg\n";
print "Max = $max\n";
print "Min = $min\n";