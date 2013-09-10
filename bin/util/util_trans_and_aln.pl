#!/usr/bin/perl -w

# -f is the fasta file of your query
# -o is the outdir
# -c is your config file
# -l is your list of accessions

#####SETUP#####

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::DB::GenBank;

my ($fasta, $outdir, $configfile, $list);
GetOptions(
	   'f|fasta=s'  => \$fasta,
	   'o|outdir=s' => \$outdir,
	   'c|config=s' => \$configfile,
	   'l|list=s'   => \$list,
	   );

#####MAIN######

`mkdir -p $outdir/pts`;
`mkdir -p $outdir/alns`;

# parse the config                                                                  
my $conf = parse_config ($configfile);

# get the fasta name
my $fasta_name;
if ($fasta =~m/\//g){
    $fasta =~m/.*\/(.*)$/;
    $fasta_name = $1;
}
else {
    $fasta_name = $fasta;
}

# get the query sequence
print STDERR "Processing query\n";
my $seqin   = Bio::SeqIO->new (-format => "fasta", 
			       -file   => $fasta);

open (T, ">$outdir/pts/$fasta_name");
# should only have one sequence
my $seqobj = $seqin->next_seq();
my $id     = $seqobj->display_id();
my $nucseq = $seqobj->seq();

my $rcobj = $seqobj->revcom;
    
# translate the sequence and check trans in each frame
my $ptobj0 = $seqobj->translate(-frame=>0);
my $ptobj1 = $seqobj->translate(-frame=>1);
my $ptobj2 = $seqobj->translate(-frame=>2);
my $rcptobj0 = $rcobj->translate(-frame=>0);
my $rcptobj1 = $rcobj->translate(-frame=>1);
my $rcptobj2 = $rcobj->translate(-frame=>2);

my $bestpt;
my $min = 10000000000000;
foreach my $pobj ($ptobj0, $ptobj1, $ptobj2, $rcptobj0, $rcptobj1, $rcptobj2){
    my $pt = $pobj->seq;
    my $ct =()= $pt =~m/\*/g;
    
    if ($ct < $min){
	$bestpt = $pt;
	$min = $ct;
    }
    else {
	next;
    }
}

# print out the best translation
print T ">$id\n$bestpt\n";

# get all the other sequences
print STDERR "Processing list\n";
my $gb = new Bio::DB::GenBank(-format => 'Fasta');

open (LIST, "$list");
while (my $line = <LIST>){
    chomp $line;

    print STDERR "Working on $line\n";

    my $seq_obj = $gb->get_Seq_by_acc($line);
    my $id = $seq_obj->id();
    my $desc = $seq_obj->desc();
    my $seq = $seq_obj->seq();

    print T ">$id\t$desc\n$seq\n";
}
close (LIST);
close (T);

# aln the translation with mafft
my $alnname = mafft_run ($fasta_name, $conf, "$outdir/pts", "$outdir/alns");

##### SUBS #####

sub mafft_run{
    my $fas  = shift;
    my $conf = shift;
    my $in   = shift;
    my $out  = shift;

    # mafft                                                                       
    my $mafft = "mafft";
    ($mafft .= " $conf->{'MAFFT'}") if ($conf->{'MAFFT'});
    $mafft .= " $in/$fas";

    `$mafft 1>$out/$fas.mafft.aln 2>$out/$fas.mafft.stderr`;
    return ("$fas.mafft.aln");
}

sub parse_config{
    my $file = shift;

    my %config;
    open (F, "$file");
    while (my $line = <F>){
	chomp $line;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }
    close (F);

    return (\%config);
}
