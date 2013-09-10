#!/usr/bin/perl -w

# -i is the alignment infile
# -f is the input format

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile, $format);
GetOptions(
	   'i|infile=s'      => \$infile,
	   'f|format=s'      => \$format,
	   );

#####MAIN#####

my $alnin = Bio::AlignIO->new(-file   => "$infile",
			      -format => "$format");

# get aln data
my $seqcounter = 0;
my $alnobj = $alnin->next_aln();
foreach my $seq ($alnobj->each_seq){
    my $id        = $seq->display_id;
    my $sequence = $seq->seq;
    
    open (W, ">$id.fa");
    print W ">$id\n$sequence\n";
    close (W);
}
