#!/usr/bin/perl -w

# -i is the directory individual alignment infiles
# -f is the input format
# -l is the list of alns you want

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($indir, $format, $outdir, $list);
GetOptions(
	   'i|indir=s'       => \$indir,
	   'f|format=s'      => \$format,
	   'l|list=s'   => \$list,
	   );

#####MAIN#####

# readin the files
opendir (D, "$indir");
my @alns = sort(readdir (D));
shift @alns;
shift @alns;
closedir (D);

my %list;
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    $list{$line} = 1;
}
close (L);

# parse the alignment
foreach my $aln (@alns){
#    my @alnname = split (/\./, $aln);
#    next unless (exists ($list{$alnname[0]}));
    next unless (exists ($list{$aln}));
    my $alnin = Bio::AlignIO->new(-file   => "$indir/$aln",
				  -format => "$format");

    # get aln data
    my $seqcounter = 0;
    my $alnobj = $alnin->next_aln();
    foreach my $seq ($alnobj->each_seq){
	$seqcounter++;
	my $id        = $seq->display_id;
	my $sequence = $seq->seq;

	(print ">$id.$aln\n$sequence\n") if ($id =~m/dmel/);
    }
}
