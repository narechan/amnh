#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use TreeSupports;
use Bio::AlignIO;
use Bio::SeqIO;

my ($help, $matrixfile, $listfile, $outdir);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'l|list=s'        => \$listfile,
    'o|outdir=s'      => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $listfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir`;

#####MAIN#####

# instantiate the object and load stuff we need
my $supportobj = TreeSupports->new;
$supportobj->load_aln    ($matrixfile);

# get the charsets
my $charsets = $supportobj->get_charsets;
my $chars    = $supportobj->get_nchar;

# store alignment information
my $partitions = {};
my $lengths    = {};
my $alnin = Bio::AlignIO->new(-file   => "$matrixfile",
			      -format => "nexus");

# store list of the partitions to keep
my $listhash = {};
open (L, "$listfile");
while (my $line = <L>){
    chomp $line;
    $listhash->{$line} = 1;
}
close (L);

# only one aln there
my $taxa = 0;
my $alnlentot = 0;
my $aln = $alnin->next_aln();
foreach my $seq ($aln->each_seq){
    my $lentot = 0;
    $taxa++;
    my $id        = $seq->display_id;
    
    foreach my $charset (sort keys %$charsets){
#	print STDERR "Storing\t$id\t$charset\n";
	next unless (exists ($listhash->{$charset}));

	my $coords = $charsets->{$charset};
	my ($start, $end) = split (/\-/, $coords);
	
	my $partition = $seq->subseq($start, $end);
	my $partlen   = length ($partition);
	$lentot += $partlen;

	$partitions->{$charset}->{$id} = $partition;
	$lengths->{$charset} = $partlen;
    }
    $alnlentot = $lentot;
}

# sort print the concatenation and charpars                                   
open (CAT, ">$outdir/matrix");
open (PRT, ">$outdir/partitions");
print CAT "#NEXUS\n";
print CAT "BEGIN DATA;\n";
print CAT "DIMENSIONS NTAX=$taxa NCHAR=$alnlentot;\n";
print CAT "FORMAT INTERLEAVE SYMBOLS=\"ACTG\" DATATYPE=DNA MISSING=? GAP=-;\n";
print CAT "MATRIX\n";
print PRT "BEGIN SETS;\n";

my $start = 1;
my $end;
foreach my $count (sort keys %$partitions){
    $end = $start - 1 + $lengths->{$count};
#    print CAT "[Partition $count length $alnlens->{$count} chars $start-$end]\n"
    ;
    print CAT "[Partition $count length $lengths->{$count} chars $start-$end]\n";
#    print PRT "CHARSET $count=$start-$end;\n";
    print PRT "CHARSET $count=$start-$end;\n";
    foreach my $sp (sort keys %{ $partitions->{$count} }){
        print CAT "$sp\t$partitions->{$count}->{$sp}\n"; 
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

