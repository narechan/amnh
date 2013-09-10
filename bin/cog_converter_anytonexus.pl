#!/usr/bin/perl -w

# -d is the directory individual alignment infiles
# -c is the cog assignment file
# -x is a string of functional categories you want
# -f is the input format
# -o is the output dir
# -t is a taxa list (in case some alignments lack certain taxa
#    and we need to insert missing blocks)

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;

my ($indir, $format, $outdir, $taxalist, $coglist, $funcstring);
GetOptions(
	   'd|indir=s'       => \$indir,
	   'f|format=s'      => \$format,
	   'o|outdir=s'      => \$outdir,
           't|taxalist=s'    => \$taxalist,
	   'c|coglist=s'      => \$coglist,
	   'x|funcstring=s'  => \$funcstring,
	   );

`mkdir -p $outdir`;

#####MAIN#####

# read in the cog data file
my %coglist;
open (C, "$coglist");
while (my $line = <C>){
    chomp $line;
    my ($part, $coggrp) = split (/\t/, $line);
    $coglist{$part} = $coggrp;
}
close (C);

# find the alignments you want
my %cogkeepers;
foreach my $cog (keys %coglist){
    if ($coglist{$cog} =~m/[$funcstring]/){
	$cogkeepers{$cog} = 1;
    }
    elsif ($coglist{$cog} eq ""){
	$cogkeepers{$cog} = 1;
    }
    else {
	next;
    }
}


# readin the files
opendir (D, "$indir");
my @alns = sort(readdir (D));
shift @alns;
shift @alns;
closedir (D);

# read in the taxa list
my @taxalist;
open (F, "$taxalist");
while (my $line = <F>){
    chomp $line;
    push (@taxalist, $line);
}
close (F);

# parse the alignment
my $alndata = {};
my $alnlens = {};
my $alnlentot = 0;
my %taxa;
foreach my $aln (@alns){
    my @alnname = split (/\./, $aln);
    next unless (exists($cogkeepers{$alnname[0]}));
    
    my $alnin = Bio::AlignIO->new(-file   => "$indir/$aln",
				  -format => "$format");

    # get aln lengths
    my $alnobj = $alnin->next_aln();
    my $alnlen = $alnobj->length;
    $alnlentot += $alnlen;
    $alnlens->{"$alnname[0]"} = $alnlen;

    # get aln data
    my %seq;
    foreach my $seq ($alnobj->each_seq){
	my $id        = $seq->display_id;
	$seq{$id} = 1;
	$taxa{$id} = 1;
	my $sequence = $seq->seq;
	$alndata->{"$alnname[0]"}->{$id} = $sequence;
    }
    
    # fill in missing data if any 
    foreach my $tax (@taxalist){
	if (exists ($seq{$tax})){
	    next;
	}
	else {
	    my $missing = '?' x $alnlen;
	    $alndata->{"$alnname[0]"}->{$tax} = $missing;
	}
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

