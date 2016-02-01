#!/usr/bin/perl -w

# -a is the input aln file
# -f is the input format
# -b is the ouput aln file
# -o is the output format

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($infile, $outfile, $informat, $outformat);
GetOptions(
	   'a|infile=s'      => \$infile,
	   'b|outfile=s'     => \$outfile,
	   'f|informat=s'    => \$informat,
	   'o|outformat=s'   => \$outformat,
	   );

my $alnin = Bio::AlignIO->new(-file   => "$infile",
			      -interleaved => 0,
			      -format => "$informat");
my $alnout = Bio::AlignIO->new(-file   => ">$outfile",
			       -interleaved => 0,
			       -format => "$outformat",
			       -idlength => 65); # hack for longer phylip names

print STDERR "$infile\n";

while (my $aln = $alnin->next_aln){
#    $aln->{'_id'} = $infile;  # hack for putting ids into stockholm files
    $alnout->write_aln($aln);
}
