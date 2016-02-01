#!/usr/bin/perl -w

# -i is the directory individual alignment infiles
# -o is the outdir

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($indir, $outdir);
GetOptions(
	   'i|indir=s'       => \$indir,
    'o|outdir=s' => \$outdir,
	   );

`mkdir -p $outdir`;

#####MAIN#####

# readin the files
opendir (D, "$indir");
my @alns = sort(readdir (D));
shift @alns;
shift @alns;
closedir (D);

# parse the alignment
my $alndata = {};
my $alnlens = {};
my $alnlentot = 0;
my %taxa;
foreach my $aln (@alns){
    my $alnin = Bio::AlignIO->new(-file   => "$indir/$aln",
				  -format => "fasta");
    
    # get aln lengths
    my $alnobj = $alnin->next_aln();
    my $alnlen = $alnobj->length;
    $alnlentot += $alnlen;
    $alnlens->{$aln} = $alnlen;
    
    # get aln data
    foreach my $seq ($alnobj->each_seq){
	my $id        = $seq->display_id;
	$taxa{$id} = 1;
	my $sequence = $seq->seq;
	$alndata->{$aln}->{$id} = $sequence;
    }
}

# sort print the concatenation and charpars                                   
my @taxa = keys %taxa;
my $taxa = @taxa;
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$taxa NCHAR=$alnlentot;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

my $start = 1;
my $end;
foreach my $count (sort keys %$alndata){
    $end = $start - 1 + $alnlens->{$count};
    print CAT "[Partition $count length $alnlens->{$count} chars $start-$end]\n";
    print PRT "CHARSET $count=$start-$end;\n";
    foreach my $sp (sort keys %{ $alndata->{$count} }){
	print CAT "$sp\t$alndata->{$count}->{$sp}\n"; 
    }
    print CAT "\n";
    $start = $end + 1;
}

print CAT ";\n";
print CAT "END;\n\n";
print PRT "END;\n";
close (CAT);
close (PRT);

`cat $outdir/matrix $outdir/partitions > $outdir/nexus`;
#`rm -rf $outdir`;


