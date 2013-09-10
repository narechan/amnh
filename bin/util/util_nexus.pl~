#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use TreeSupports;
use Bio::AlignIO;
use Bio::SeqIO;

my ($help, $matrixfile, $outdir);
GetOptions(
    'h|help'          => \$help,
    'm|matrix=s'      => \$matrixfile,
    'o|outdir=s'      => \$outdir,
	   ) or pod2usage;
pod2usage if $help;

for my $option ($matrixfile, $outdir){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

`mkdir -p $outdir/fastas`;

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

# only one aln there
my $aln = $alnin->next_aln();
foreach my $seq ($aln->each_seq){
    my $id        = $seq->display_id;
    
    foreach my $charset (sort keys %$charsets){
#	print STDERR "Storing\t$id\t$charset\n";

	my $coords = $charsets->{$charset};
	my ($start, $end) = split (/\-/, $coords);
	
	my $partition = $seq->subseq($start, $end);
	my $partlen   = length ($partition);
	
	$partitions->{$charset}->{$id} = $partition;
	$lengths->{$charset} = $partlen;
    }
}

# print out fasta file alignments of each charset
foreach my $cset (sort keys %$partitions){
    open (FH, ">$outdir/fastas/$cset.fa");

    foreach my $id (sort keys %{ $partitions->{$cset} }){
        my $seq = $partitions->{$cset}->{$id};
	my $len = length ($seq);
#	my $matchstring = '\\' . '?' . "{$len}";
	my $matchstring = '\\' . '-' . "{$len}";
	
	next if ($seq =~m/$matchstring/);

#	my $newid;
#	if ($id =~m/^2/){
#	    $newid = "Sa_" . $id;
#	}
#	else{
#	    $newid = $id;
#	}
	
#	print FH ">$newid\n$seq\n";
	print FH ">$id\n$seq\n";
    }
    close (FH);
}

# cycle through the characters and charsets and count missing
my %freq;
foreach my $cset (sort keys %$partitions){
    my $len    = 0;
    $freq{$cset} = 0;

    foreach my $id (sort keys %{ $partitions->{$cset} }){
	my $seq = $partitions->{$cset}->{$id};
	$len = length ($seq);

	# full taxa per partition
	my $matchstring = '\\' . '-' . "{$len}"; # a bit hokey!                   
	($freq{$cset}++) unless ($seq =~m/$matchstring/);
    }
}

open (F, ">$outdir/freqs.stats");
foreach my $freq (sort keys %freq){
    print F "$freq\t$freq{$freq}\n";
}
close (F);

