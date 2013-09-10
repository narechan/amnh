#!/usr/bin/perl -w

# -g is the directory of input gff3 annotations
# -i is the input nexus aln
# -l is the lookup file of OG to accesssions
# -o is the outdir

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SeqIO;
use TreeSupports;

my ($inalign, $gff3dir, $lookupfile, $outdir);
GetOptions(
	   'i|inalign=s'    => \$inalign,
	   'g|gff3dir=s'    => \$gff3dir,
	   'l|lookup=s'       => \$lookupfile,
	   'o|outdir=s'       => \$outdir,
	   );

#####MAIN######

`mkdir -p $outdir`;

# get the nuc accs for each partition using the lookup file
my $lookup = {};
open (L, "$lookupfile");
while (my $line = <L>){
    chomp $line;
    my ($og, $accs) = split (/\t/, $line);
    my @accs = split (/\,/, $accs);
    $lookup->{$og} = [@accs];
}
close (L);

# get all available annotations from the gff3 files
opendir (D, "$gff3dir");
my @gff3s = sort (readdir (D));
shift @gff3s;
shift @gff3s;
closedir (D);

my $annotations = {};
foreach my $gff3 (@gff3s){
    my ($sp, $crap) = split (/\./, $gff3);

    open (F, "$gff3dir/$gff3");
    while (my $line = <F>){
	next if ($line =~m/^\#/);
	chomp $line;
	
	my @line = split (/\t/, $line);
	my @props = split (/\;/, $line[8]);
	
	my $id;
	my $annot;
	foreach my $prop (@props){
	    my ($key, $value) = split (/\=/, $prop);
	    
	    if ($key eq "ID"){
		$value =~s/fig\|//g;
		$value =~s/\.//g;
		$id = "Sa" . $sp . "#" . $value;
	    }
	    elsif ($key eq "Name"){
		$annot = $value;
	    }
	    else{
		next;
	    }
	}
	
#	print STDERR "$id\t$annot\n";
	
	$annotations->{$id} = $annot;
    }
}

# update the nexus file with group members and annotations 
# as comments
my $supportobj = TreeSupports->new;
$supportobj->load_aln    ($inalign);
my $charsets = $supportobj->get_charsets;
my $chars    = $supportobj->get_nchar;

my $partitions = {};
my $lengths    = {};
my $alnin = Bio::AlignIO->new(-file   => "$inalign",
                              -format => "nexus");

my $aln = $alnin->next_aln();

my $alndata = {};
my $alnlens = {};
my $alnlentot = 0;
my @taxa;
foreach my $seq ($aln->each_seq){
    push (@taxa, $seq);
    my $id        = $seq->display_id;

    my $lentot = 0;
    foreach my $charset (sort keys %$charsets){

        my $coords = $charsets->{$charset};
        my ($start, $end) = split (/\-/, $coords);

        my $partition = $seq->subseq($start, $end);
	my $partlen   = length ($partition);
	$lentot += $partlen;

	$alnlens->{$charset} = $partlen;
	$alndata->{$charset}->{$id} = $partition;
    }
    $alnlentot = $lentot;
}

# sort print the concatenation and charpars                                       
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
    my @accs = @{$lookup->{$count}};
    my @annots;
    foreach my $acc (@accs){
	my $annotst;
	if ($annotations->{$acc}){
	    $annotst = $acc . "=" . $annotations->{$acc};
	}
	else {
	    $annotst = $acc . "=" . "NULL";
	}
	push (@annots, $annotst);
    }
    my $astring = join ";", @annots;
    print CAT "[$astring]\n";
    
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
